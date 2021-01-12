close all
clc
clear all
% chirp variables
SF = 8;
% step_size = 1e3;
N = 2^SF;
BW = 125e3;%(2^SF) * step_size;
Fs = 125e3;%sampling_freq_factor*BW;
sampling_freq_factor = Fs/BW;
upsampling_factor = Fs/BW;
Ts = 1/Fs;
t = 0:Ts:1;
SNR = 30;

Freq_off = 0;
% Chirplet Transform variables
sigma = 0.5;%0.94;%0.31623;%0.288675;
fLevel = N;
WinLen = N;
alpha = -(BW^2)/(2^SF);%128e6; % in Hz/s¿


% generating LORA pkts
num_preamble = 8;
num_sync = 2;
preamble_sym = 1;
num_data_sym = 60;
num_DC = 2.25;
pkt_len = num_preamble + num_DC + num_data_sym + num_sync;
num_samples = pkt_len * N * (Fs/BW);
num_TX = 5;%j;%20;


% UC = sym_to_data(1,N);
UC = sym_to_data_ang([1 1 1],N);
DC = conj(UC);
DC_SFD = DC(1:2.25*(N*sampling_freq_factor));
DC_pre = conj(sym_to_data_ang( preamble_sym  .*  ones(1,num_preamble),N));

sync = sym_to_data_ang([9 17],N);
start_frame = [sym_to_data_ang(ones(1,num_preamble),N) sync DC_SFD];

% sym(1,:) = [1:1+num_data_sym-1];
% sym(2,:) = [20:20+num_data_sym-1];
a = 0.3;
b = 1;
for i = 1:num_TX
%     rng(i);
    r = (b-a).*rand(1,1) + a
%     rng(1);
    sym(i,:) = [randi([1 N],1,num_data_sym)];%[1 6 11 16 21 26 31 36 41 46];
    Data(i,:) = [start_frame sym_to_data_ang(sym(i,:),N)];
%     c_data = c_data + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating random offset for each transmitter and adding the signals
% rng(20);
offset = randi(N* (pkt_len),num_TX,1);%[20]%; 300+N+7];%[300; 330];%randi(N* (pkt_len),num_TX,1);% [100;71;1270];% [300;304;308;312;316];%    %

% if( k == 3)
%         keyboard
% end

% offset = sortrows(offset);
% if(num_TX > 1)
%     b = nchoosek(offset(:),2);
%     while(min(mod(abs(b(:,1) - b(:,2)),N)) < 10 || max(mod(abs(b(:,1) - b(:,2)),N)) > 246)
%         offset = randi(N* (pkt_len),num_TX,1);%
%         offset = sortrows(offset);
%         b = nchoosek(offset(:),2);
% %         disp('Stuck - exit please')
%     end
%     min(mod(abs(b(:,1) - b(:,2)),N))
%     max(mod(abs(b(:,1) - b(:,2)),N))
% end

% for i = 1:num_TX
%     offset(i,:) = 300 + (i-1)*4;
% end

o = sort(offset)+1 + ((num_preamble + num_sync)*N);
sort_offset = sort(offset);
o(:,2) = o(:,1) + N;

% min(abs(b(:,1) - b(:,2)))
% o

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0)
    Data(1,:) = Data(1,:) .* exp((j*2*pi*(Freq_off/Fs)).*(1:length(Data(1,:))));
end
if(0)
    Data(1,:) = Data(1,:) .* exp((j*2*pi*(2/128)).*(1:length(Data(1,:))));
end

temp = zeros(num_TX,num_samples + max(offset(1:num_TX)));
for i = 1:num_TX
    temp(i,offset(i)+1 : offset(i)+num_samples) = Data(i,:);%     zeros(1,2*N)
end

Rx_Buffer = sum(temp,1);

Rx_Buffer = awgn(Rx_Buffer,SNR,'measured');
Rx_Buffer = [Rx_Buffer zeros(1,2*N)];

off_sym = [offset sym];
off_sym = sortrows(off_sym);
sym = off_sym(:,2:end);
%%      DC correlations
ind = DC_location_correlation(Rx_Buffer(1:upsampling_factor:end),N,DC(1:N),16,0.2);

temp = [];
indices = [zeros(1,floor(num_DC)); ind];

for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) > 3 )
%         temp = [temp; indices(i,:)];
%     end
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 3 )
            temp = [temp; indices(i,:)];
        end
    end
end
DC_ind = temp

%%      UC correlations
Rx_Buffer = Rx_Buffer(1:floor(length(Rx_Buffer)/upsampling_factor)*upsampling_factor);
for i = 1:upsampling_factor
    Rx_Buff_dnsamp(i,:) = Rx_Buffer(i:upsampling_factor:end);
end

% [Data_freq_off Upchirp_ind] = UC_location_corr(Rx_Buffer(1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind);
[Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer(1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind,16,0.15);
temp = [];
indices = [zeros(1,num_preamble); Upchirp_ind];

for i = 2:size(indices,1)
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 5 )
            temp = [temp; indices(i,:)];
        end
    end
end
Upchirp_ind = temp

[Data_freq_off, Peak] = dnsamp_buff(Rx_Buff_dnsamp,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre);
[Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre,Fs,Peak);


temp = [];
temp_data = [];
temp_peaks = [];
indices = [zeros(1,num_preamble); Preamble_ind];
Data = [zeros(1,size(Data_out,2)); Data_out];
peaks = [zeros(1,size(Peak_amp,2)); Peak_amp];
clear Peak_amp
for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) <= 3 )
%         
%     else
%         temp = [temp; indices(i,:)];
%         temp_data = [temp_data; Data(i,:)];
%     end

    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
        temp_data = [temp_data; Data(i,:)];
        temp_peaks = [temp_peaks; peaks(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 5 )
            temp = [temp; indices(i,:)];
            temp_data = [temp_data; Data(i,:)];
            temp_peaks = [temp_peaks; peaks(i,:)];
        end
    end
end
Pream_ind = temp
Data_out = temp_data;
Peak_amp = temp_peaks;
%%  Data Demodulations
demod_sym = [];
% num_data_sym = 43;%63;%39;
for j =1:size(Pream_ind,1)
    j
    [demod_sym(j,:) sym_peak(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
    demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j),N);
end
%%
% load('SF8/SF8_symbols.mat');
% num_TX = size(Pream_ind,1);
% sym = sym(1:num_data_sym);
% sym = repmat(sym,num_TX,1);
ser = ((num_data_sym * num_TX) - sum(sum(demod_sym == sym)))/(num_data_sym * num_TX);
pkt_ser = num_data_sym - sum(demod_sym == sym,2)
ser*(num_TX*num_data_sym)
ser

%%
i = 5;
ind = find(demod_sym(i,:) ~= sym(i,:))
sym(i,ind)
demod_sym(i,ind)
%% Packet interference
interf_map = zeros(2*size(Pream_ind,1) + 2,length(Rx_Buff_dnsamp));
c = 1;
for i = 2 : 2 : 2*size(Pream_ind,1)
    interf_map(i, Pream_ind(c,1) : Pream_ind(c,1) + (num_preamble+num_sync) * N - 1 ) = (3/4)*ones(1,(num_preamble+num_sync) * N);
    interf_map(i, Pream_ind(c,1) + (num_preamble+num_sync) * N : Pream_ind(c,1) + (num_preamble+num_sync+num_DC) * N - 1 ) = (2/4)*ones(1,(num_DC) * N);
    interf_map(i, Pream_ind(c,1) + (num_preamble+num_sync+num_DC) * N : Pream_ind(c,1) + (num_preamble+num_sync+num_DC+num_data_sym) * N - 1 ) = ones(1,(num_data_sym) * N);
    c = c+1;
end
spec_plot(interf_map,size(interf_map,1)-1,0,0,1);
% plot(interf_map(1,:))
% hold on 
% plot(interf_map(2,:))

