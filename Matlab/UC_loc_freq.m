function [Data_freq_off in ind] = UC_location(Data,N,UC,o,num_preamble,num_sync,DC,BW,SF,sigma,Fs)
%UC_LOCATION Summary of this function goes here
%   Detailed explanation goes here
dis = 0;
fLevel = N;
WinLen = N;
alpha = 0;
% sigma = 0.1;

Upchirp_ind = [];
edge_range_threshold = 4;

c = 1;
if(size(o,1) == 0)
    return;
end
for i = 1:size(o,1)
    if(o(i,1) - ((num_preamble+num_sync)*N) < 1)
        continue;
    end
    pot_pream_ind(c,:) = o(i,1) - ((num_preamble + num_sync)*N) - 1 : N : o(i,1)- ((num_sync)*N);
    c = c+1;
end

count = 0;
max_fft = [];
max_fft_ind = [];

for k = 1:size(pot_pream_ind,1)
        close all
        v = 15;
        for i = -v:v

            temp_wind(i+v+1) = sum(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC)...
            ./ sqrt(sum( Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i) .* conj(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i)) ) .* ...
            sum( DC .* conj(DC)));
        end
        if(dis == 1)
            figure
            plot([-v:v],abs(temp_wind))
        end
        
        [~, b] = max(abs(temp_wind));
        b = b - v - 1;
        b = floor((b)/2);
        data_wind = Data(pot_pream_ind(k,1) + 1 + b:pot_pream_ind(k,num_preamble) + N + b); %pot_pream_ind(k,num_preamble) + N);
        
        
        freq_off = [];
        for i = 1:num_preamble
            data_fft(i,:) = abs(fft(data_wind((i-1)*N + 1:i*N).*DC(1:N),128*N));
            [~,c(i)] = max(data_fft(i,:));
            if(c(i) > (128*N)/2)
                freq_off = [freq_off ((N*128) - c(i))/128];
            else
                freq_off = [freq_off -(c(i) - 1)/128];
            end
        end
        
        freq_off = sum( freq_off(2:7) ) / (num_preamble - 2);
        freq_off = freq_off - floor(freq_off);
%         for i = 1:num_preamble
% 
%             data_fft(i,:) = abs(fft(data_wind((i-1)*N + 1:i*N).*DC(1:N),128*N));
% %             data_fft_org(i,:) = abs(fft(data_wind_org((i-1)*N + 1:i*N).*DC(1:N),500*N));
%         end
%         [~,c] = max(data_fft(2,:));
%         freq_off = ((N*128) - c)/128
        if(dis == 1)
            keyboard
        end
        Data_freq_off = Data .* exp((j*2*pi*(freq_off./N)).*(1:size(Data,2)));
        pot_pream_ind = pot_pream_ind + b;
        ind(k,:) = o(k,:) + b;
        
        for i = -v:v
%             temp_wind(i+20+1,:) = abs(fft(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC));
%             temp_wind(i+20+1) = sum(Data(pot_pream_ind(k,1) + i + 1:pot_pream_ind(k,2) + i).*DC(1:N));
%             temp_wind(i+1) = sum(Data(pot_pream_ind(k,1) + i + 1:pot_pream_ind(k,2) + i).*DC(1:N))...
%             / sqrt(sum( Data(pot_pream_ind(k,1) + i + 1:pot_pream_ind(k,2) + i) .* conj(Data(pot_pream_ind(k,1) + i + 1:pot_pream_ind(k,2) + i)) ) * ...
%             sum( DC(1:N) .* conj(DC(1:N))));
            temp_wind_new(i+v+1) = sum(Data_freq_off(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC)...
            ./ sqrt(sum( Data_freq_off(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i) .* conj(Data_freq_off(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i)) ) .* ...
            sum( DC .* conj(DC)));
%             temp_wind(i+1) = sum(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC);
%             max_fft = [max_fft max(temp_wind(i+20+1,:))];
        end
%         subplot(313)
        if(dis == 1)
            figure
            plot([-v:v],abs(temp_wind_new))
        end
        data_wind_freq = Data(pot_pream_ind(k,4) + 1:pot_pream_ind(k,4) + N);
        data_fft = abs(fft(data_wind_freq.*DC(1:N),N));
%         [~,pre] = max(data_fft);
%         pre
        
%         in = abs(temp_wind_new(v + 1));
        in = max(abs(temp_wind_new));
        if(dis == 1)
            keyboard
        end
        Upchirp_ind = [Upchirp_ind; pot_pream_ind(k,:) + 1];
        
end

end