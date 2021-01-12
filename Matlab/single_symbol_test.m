%%      First Section
close all;
clc
clear all;

SF = 7;
% preamble_sym = 1;
% num_preamble = 1;
Freq_off = 1700;
% overlap = 0;

step_size = 1000;
BW = (2^SF) * step_size;
sampling_freq_factor = 1;
Fs = sampling_freq_factor*BW;
Ts = 1/Fs;
chirp_rate = BW;
frame = ceil(Fs/chirp_rate);
N = frame * (2^SF);
SNR = 40;

% Chirplet Transform variables
% sigma = 0.05;%0.94;%0.94;%0.94;%0.31623;%0.288675;
edge_range_threshold = 10;
fLevel = N;
WinLen = N;
alpha = 10;%(BW^2)/(2^SF);%128e6; % in Hz/s¿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbols1 = [24];                                             % Symbols of TX1
symbols2 = [34 65];                                           % Symbols of TX2
symbols3 = [24 49];
symbols4 = [115 102];

Data1 = sym_to_data_ang(symbols1,N);           %   TX1
Data2 = sym_to_data_ang(symbols2,N);
% Data3 = sym_to_data_ang(symbols3,N);
% Data4 = sym_to_data_ang(symbols4,N);
DC = conj(sym_to_data_ang(1,N));

i = 20;
i1 = 30;
i2 = 40;

window = Data1 + Data2(N - i:(N-i) +N - 1); % Data1(N+1 : 2*N )%[zeros(1,10) Data1 zeros(1,10)];%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
    Rx_Buffer = window .* exp((j*2*pi*(Freq_off/Fs)).*(1:length(window)));
end
% window = Data1 + Data2(N - i:(N-i) +N - 1) + Data3(N - i1:(N-i1) +N - 1);% Data1(N+1 : 2*N )
Rx_Buffer = awgn(window,SNR,'measured');


c = 1;
wind = [];
out = [];
for sigma = [5 0.05]
%         close all
    [temp,~,~] = Chirplet_Transform(Rx_Buffer.*DC,fLevel,WinLen,Fs,0,sigma);
     wind(c,:,:) = abs(temp);
%         wind(c,:,:) = (abs(temp).^2)./max(max(abs(temp).^2));
%         wind_edge(c,:,:) = edge(abs(temp));
%     if(dis)
%         spec_plot(abs(temp),N,0,1,0,0)
%     end
%         spec_plot(temp,N,0)
    c = c + 1;
end
for i = 1:N+1
    for j = 1:size(wind,3)
        absSpec = abs(wind(:,i,j));
        index = find(min(absSpec) == absSpec);
        out(i,j) = wind(index(1),i,j);
    end
end
spec_plot(out,N,0,0,1);
Rx_Buffer = [zeros(1,64) Rx_Buffer];
DC = conj(sym_to_data_ang(65,N));
stft(Rx_Buffer,N,DC,0);
% [temp,~,~] = Chirplet_Transform(Rx_Buffer.*DC,fLevel,WinLen,Fs,alpha,0.2);
% spec_plot(abs(temp).^2,N,0,0,0)

% %Cross-Corr
% for i = 1:length(Rx_Buffer) - N
% %             temp_wind(i+20+1,:) = abs(fft(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC));
% %             temp_wind(i+20+1) = sum(Data(pot_pream_ind(k,1) + i + 1:pot_pream_ind(k,2) + i).*DC(1:N));
%     temp_wind(i) = sum(Rx_Buffer(i:N+i-1).*DC)...
%     / sqrt(sum( Rx_Buffer(i:N+i-1) .* conj(Rx_Buffer(i:N+i-1)) ) * ...
%     sum( DC .* conj(DC)));
% %             temp_wind(i+20+1) = sum(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC);
% %             temp_wind(i+1) = sum(Data(pot_pream_ind(k,1) + 1 + i:pot_pream_ind(k,num_preamble) + N + i).*DC);
% %             max_fft = [max_fft max(temp_wind(i+20+1,:))];
% end
% 
% plot(abs(temp_wind))
% 
% %Acorr
