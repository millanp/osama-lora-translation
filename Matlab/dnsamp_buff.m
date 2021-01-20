function [Data_buff peak_amp Up_ind] = dnsamp_buff(Data_stack,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC)
%DNSAMP_BUFF Summary of this function goes here
%   Detailed explanation goes here

Up_ind = [];
peak_amp = [];
Data_buff = [];
n_pnt = 16;
for k = 1:size(Upchirp_ind,1)
    k
    if(Upchirp_ind(k,1) - N <= 0)
        continue;
    end
    close all
    in = [];
    for m = 1:size(Data_stack,1)
        m
        %         pream_fft = abs(fft(Data_stack(m,Upchirp_ind(k,1) : Upchirp_ind(k,2) - 1) .* DC(1:N)));
        %         [~,pream_bin] = max(pream_fft);
        %         dnchirp_fft = abs(fft(Data_stack(m,Upchirp_ind(k,1) + (num_preamble+num_sync)*N : Upchirp_ind(k,1) + ((num_preamble+num_sync)*N) + N - 1) .* conj(DC(1:N))));
        %         [~,DC_bin] = max(dnchirp_fft);
        %         figure
        %         plot(pream_fft)
        %         figure
        %         plot(dnchirp_fft)
        %         keyboard
        data_wind = [];
        data_fft = [];
        freq_off = [];
        ind_temp = [1:5*n_pnt (N*n_pnt)-(4*n_pnt):(N*n_pnt)];
        for j = 1:num_preamble
            data_wind = Data_stack(m,Upchirp_ind(k,1) : Upchirp_ind(k,1) + (num_preamble*N) -1);
            data_wind_slice = data_wind((j-1)*N + 1:j*N);
            first_fft_arg = data_wind((j-1)*N + 1:j*N) .* DC(1:N);
            data_fft(j,:) = abs(fft(data_wind((j-1)*N + 1:j*N) .* DC(1:N),n_pnt*N));
            
            %                 plot(data_fft(j,:))
            
            [~,c(j)] = max(data_fft(j,ind_temp));
            c(j) = ind_temp(c(j));
            if(c(j) > (n_pnt*N)/2)
                freq_off = [freq_off ( (N*n_pnt) - c(j) ) / n_pnt];
            else
                freq_off = [freq_off -1*( c(j) - 1 ) / n_pnt];
            end
            
        end
        freq_off = sum( freq_off(2:7) ) / (num_preamble - 2);
        Data_freq_off(m,:) = Data_stack(m,:) .* exp( (1i*2*pi*(freq_off./N)) .* (1:length(Data_stack(m,:))) );
        
        clear data_wind data_fft ind_temp
        ind_temp = [1:5 (N-4):N];
        a = [];
        for j = 1:num_preamble
            data_wind = Data_freq_off(m,Upchirp_ind(k,1) : Upchirp_ind(k,1) + (num_preamble*N) -1);
            data_fft(j,:) = abs(fft(data_wind((j-1)*N + 1:j*N) .* DC(1:N),N));
            %                 plot(data_fft(j,:))
            [a(j),c(j)] = max(data_fft(j,ind_temp));
            c(j) = ind_temp(c(j));
        end
        peak_stats(k,m,1) = mean(a);
        %             if(std(a) < 0.1)
        %                 peak_stats(k,m,2) = 0.1;
        %             elseif(std(a) > 0.3)
        %                 peak_stats(k,m,2) = 0.3;
        %             else
        peak_stats(k,m,2) = var(a);%(max(a) - mean(a)) + (mean(a) - min(a));
        peak_stats(k,m,3) = std(a);
        useless1 = peak_stats(k,m,:);
        %             end
        
        
        %         v = 15;
        %         temp_wind = [];
        %         for i = -v:v
        %             temp_wind(i + v + 1) = sum(Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1) .* DC)...
        %             ./ sqrt(sum( Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1) .* conj(Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1)) ) .* ...
        %             sum( DC .* conj(DC)));
        %         end
        %         figure
        %         plot([-v:v],abs(temp_wind))
        %         keyboard
        specarg1 = Data_freq_off(m,Upchirp_ind(k,1) - N:Upchirp_ind(k,end) + N - 1 - N);
        Spec = stft(Data_freq_off(m,Upchirp_ind(k,1) - N:Upchirp_ind(k,end) + N - 1 - N),N,DC(1:N),0,0);
        temp = [];
        freq_track_qual = [];
        pream_peak_ind = [];
        adj_ind = [];
        row_ind = [N-5:N 1:6];
        count = 1;
        for i = row_ind
            temp(count) = sum(abs(Spec(i,:)));
            count = count + 1;
        end
        [~,ind] = max(temp);
        pream_peak_ind = row_ind(ind);
        adj_ind = [mod(pream_peak_ind-1,N) mod(pream_peak_ind+1,N)];
        if(sum(adj_ind == 0) == 1)
            adj_ind(find(adj_ind == 0)) = N;
        end
        freq_track_qual = ( sum(abs(Spec(pream_peak_ind,:))) - sum(abs(Spec(adj_ind(1),:))) ) + ( sum(abs(Spec(pream_peak_ind,:))) - sum(abs(Spec(adj_ind(2),:))) );
        %             keyboard
        in = [in freq_track_qual];
        
        
        %         in = [in max(abs(temp_wind))];
        
    end
    [~,b] = max(in);
    %     Data_buff(k,:) = Data_freq_off(b,:);
    Data_buff = [Data_buff; Data_freq_off(b,:)];
    %     peak_amp(k,:) = reshape(peak_stats(k,b,:),1,[]);
    peak_amp = [peak_amp; reshape(peak_stats(k,b,:),1,[])];
    Up_ind = [Up_ind; Upchirp_ind(k,:)];
    
end

end
