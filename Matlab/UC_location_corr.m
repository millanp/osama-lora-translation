function [Data_freq_off Upchirp_ind] = UC_location_corr(Data,N,num_preamble,num_sync,num_DC,num_data_sym,DC,BW,SF,Fs,DC_ind)
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

for i = 1:length(Data) - length(DC)
    temp_wind(i+1) = sum(Data(i + 1 : i + N).*DC(1:N))...
    / sqrt(sum( Data(i + 1: i + N) .* conj(Data(i + 1 : i + N)) ) * ...
    sum( DC(1:N) .* conj(DC(1:N))));
end
plot(abs(temp_wind));
% keyboard

n_samp_array = [];
peak_ind_prev = [];
for i = 0:floor(length(temp_wind)/N)-1
    wind = abs(temp_wind(i*N + 1 : (i+1) * N));
    peak_ind_curr = get_4_max(wind,0.2,16);
    
    if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)
        
        for j = 1:length(peak_ind_curr)
            for k = 1:length(peak_ind_prev)
                
                if(abs(peak_ind_curr(j) == peak_ind_prev(k)))
%                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                    n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N)];
                end
                
            end
        end
        
    end
    
    peak_ind_prev = peak_ind_curr;
end

for i = 1:length(n_samp_array)
    c = 0;
    ind_arr = n_samp_array(i) + N : N : n_samp_array(i) + N + ((num_preamble-2)*N);
    
    for j = 1:length(ind_arr)
        c = c + sum( n_samp_array == ind_arr(j) );
    end

    if( c >= 6 )
        Upchirp_ind = [Upchirp_ind; [n_samp_array(i) ind_arr]];
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

Data_freq_off = 0;

end