function [Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Data,N,num_preamble,num_sync,num_DC,num_data_sym,DC,BW,SF,Fs,DC_ind,pnts_threshold,corr_threshold)
%UC_LOCATION Summary of this function goes here
%   Detailed explanation goes here

if(size(DC_ind,1) == 0)
    return;
end
pot_pream_ind = [];
c = 1;
for i = 1:size(DC_ind,1)
    if(DC_ind(i,1) - ((num_preamble+num_sync)*N) < 1)
        continue;
    end
    pot_pream_ind(c,:) = DC_ind(i,1) - ((num_preamble + num_sync)*N) : N : DC_ind(i,1)- ((num_sync)*N);
    c = c+1;
end

Upchirp_ind = [];

temp_wind = [];
for j = 1:size(pot_pream_ind,1)
    if(pot_pream_ind(j,1) - N <= 0)
        continue;
    end
    Data_buffer = [];
    Data_buffer = Data(pot_pream_ind(j,1) - N : pot_pream_ind(j,end)-1 + N);
    temp = [];
    for i = 1:length(Data_buffer) - length(DC)
        temp(i+1) = sum(Data_buffer(i + 1 : i + N).*DC(1:N))...
        / sqrt(sum( Data_buffer(i + 1: i + N) .* conj(Data_buffer(i + 1 : i + N)) ) * ...
        sum( DC(1:N) .* conj(DC(1:N))));
    end
    temp_wind(j,:) = temp;
end
% keyboard
% figure
% plot(abs(temp_wind(1,:)));
% plot(abs(temp_wind(2,:)));

array_stack = {};
for m = 1:size(temp_wind,1)
    
    n_samp_array = [];
    peak_ind_prev = [];
    for i = 0:floor(length(temp_wind)/N)-1
        
        wind = abs(temp_wind(m,i*N + 1 : (i+1) * N));
        peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold);

        if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)

            for j = 1:length(peak_ind_curr)
                for k = 1:length(peak_ind_prev)

                    if(abs(peak_ind_curr(j) == peak_ind_prev(k)))
    %                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                        n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N)+(pot_pream_ind(m,1)-N-1)];
                    end

                end
            end

        end

        peak_ind_prev = peak_ind_curr;
    end
    array_stack{m} = n_samp_array;

end

for m = 1:length(array_stack)
    n_samp_array = [];
    n_samp_array = cell2mat(array_stack(m));
    
    for i = 1:length(n_samp_array)
        c = 0;
        ind_arr = n_samp_array(i) + N : N : n_samp_array(i) + N + ((num_preamble-2)*N);

        for j = 1:length(ind_arr)
            c = c + sum( n_samp_array == ind_arr(j) );
        end
        
        if( c >= 6 )
            if(length(Upchirp_ind) ~= 0)
                if(sum(n_samp_array(i) == Upchirp_ind(:,1)) ~= 1)
                    Upchirp_ind = [Upchirp_ind; [n_samp_array(i) ind_arr]];
                else
                    
                end
            else
                Upchirp_ind = [Upchirp_ind; [n_samp_array(i) ind_arr]];
            end
        end
    end
    
end
% c = 0;
% for i = 1:size(Upchirp_ind,1)
%     data_wind = [];
%     data_fft = [];
%     freq_off = [];
%     for j = 1:num_preamble
%         data_wind = Data(Upchirp_ind(i,1) : Upchirp_ind(i,1) + (num_preamble*N) -1);
%         data_fft(j,:) = abs(fft(data_wind((j-1)*N + 1:j*N) .* DC(1:N),128*N));
%         plot(data_fft(j,:))
%         [~,c(j)] = max(data_fft(j,:));
%         if(c(j) > (128*N)/2)
%             freq_off = [freq_off ( (N*128) - c(j) ) / 128];
%         else
%             freq_off = [freq_off -1*( c(j) - 1 ) / 128];
%         end
%     end
%     freq_off = sum( freq_off(2:7) ) / (num_preamble - 2);
%     Data_freq_off(i,:) = Data .* exp( (1i*2*pi*(freq_off./N)) .* (1:length(Data)) );
% end
temp = [];
indices = [zeros(1,num_preamble); Upchirp_ind];

for i = 2:size(indices,1)
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        x = indices(i)
        y = temp(:,1)
        thing = min(abs(indices(i) - temp(:,1)))
        if( min(abs(indices(i) - temp(:,1))) > 5 )
            temp = [temp; indices(i,:)];
        end
    end
end
Upchirp_ind = temp;

Data_freq_off = 0;

end