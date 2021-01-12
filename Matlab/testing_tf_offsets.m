close all;
% clear all;
clc
clearvars -except cursor_info cursor_info1

% chirp variables
SF = 8;
N = 2^SF;
BW = 250e3;
Fs = 2e6;
upsampling_factor = Fs/BW;
Ts = 1/Fs;

% SNR = 10;

% Chirplet Transform variables
sigma = 0.5;%0.94;%0.31623;%0.288675;
fLevel = N;
WinLen = N;
alpha = -(BW^2)/(2^SF);%128e6; % in Hz/s¿


% generating LORA pkts
num_preamble = 8;
num_sync = 2;
preamble_sym = 1;
num_data_sym = 15;
num_DC = 2.25;
pkt_len = num_preamble + num_DC + num_data_sym + 2;% + num_sync
num_samples = pkt_len * N;
% num_TX = 1;%j;%20;

DC = conj(sym_to_data_ang([1],N));

UC_SFD = sym_to_data_upsampled([1 1 1],N,Fs,BW);
UC_SFD = UC_SFD(1:2.25*8*N);
DC_upsamp = conj(sym_to_data_upsampled([1],N,Fs,BW));

DC_pre = conj( sym_to_data_ang(preamble_sym  .*  ones(1,num_preamble),N) );
DC_pre_upsamp = conj( sym_to_data_upsampled(preamble_sym  .*  ones(1,num_preamble),N,Fs,BW) );
sync = sym_to_data_ang([9 17],N);
start_frame = [sym_to_data_ang(ones(1,num_preamble),N) sync conj(UC_SFD)];

% fi_1 = fopen('adafruit_SF7');
% fi_1 = fopen('SF8/5_tx');
% fi_1 = fopen('edge_offset/1tx_SF8_imp');
fi_1 = fopen('C:\Osama\Matlab_lt\indoor_15nodes_30dBspd_out\15tx_1lm_1');
% fi_1 = fopen('C:\temp\USRP_data\15tx_attendance');
% fi_1 = fopen('5_tx_SF8');
% fi_1 = fopen('5_tx_SF8');
% fi_1 = fopen('collisions_2');
% fi_1 = fopen('collisions_usrp.dat');
x_inter_1 = fread(fi_1, 'float32');
fclose(fi_1);
% if data is complex
x_1 = x_inter_1(1:2:end) + 1i*x_inter_1(2:2:end);

x_1 = x_1.';
t = [0:length(x_1)-1]/Fs;

x_1 = x_1(1:floor(length(x_1)/upsampling_factor)*upsampling_factor);

x_1_dnsamp = x_1(1:upsampling_factor:end);

% subplot(211)
% plot(real(x_1))
% subplot(212)
plot(real(x_1_dnsamp))
hold on 
% grid minor
plot(0.01 *ones(1,length(x_1_dnsamp)));
%% extract Packets
RSSI_threshold = 1e-3;
pnts = find(x_1 > RSSI_threshold);
pkt_bndrs = [];
for i = 1:length(pnts)-1
    if(pnts(i+1) - pnts(i) > 2e3)
        pkt_bndrs = [pkt_bndrs (pnts(i) + (N*50)) (pnts(i+1) - (N*50))];
    end
end
pkt_bndrs = [(pnts(1) - (N*50)) pkt_bndrs (pnts(end) +  (50*N))];
pkts = reshape(pkt_bndrs,2,[])';
disp(['found ' num2str(size(pkts,1)) ' Packets']);
%%
Packet_number = input('Which Packet you want to proceed with : ');
% ST_LORA = pkts(:,1);
% ED_LORA = pkts(:,2);
ST_LORA = cursor_info.DataIndex*upsampling_factor;%2.18e5;%6e5;%2.385e5;%5.58e5;%1.346e6;%4.23e5;%2.38e5;%         3.43e6;%
ED_LORA = cursor_info1.DataIndex*upsampling_factor;%2.25e5;%6.16e5;%2.475e5;%5.68e5;%1.353e6;%4.3e5;%2.46e5;%          3.48e6;%

% adafruit
% ST_LORA = [3.1e4 2.9e5 5.492e5 8.08e5 1.067e6];
% ED_LORA = [3.95e4 2.99e5 5.575e5 8.162e5 1.0752e6];

Rx_Buffer = x_1(ST_LORA(Packet_number):ED_LORA(Packet_number));

% plot(real(Rx_Buffer))
plot(real(Rx_Buffer(1,1:8:end)))
%%  Novel DC Location
% i = 12159;
DC_ind = novel_DC_location(Rx_Buffer(1:8:end),N,DC,num_DC,10);
DC_ind = sortrows(DC_ind)

%%      UC correlations
Rx_Buffer = Rx_Buffer(1:floor(length(Rx_Buffer)/upsampling_factor)*upsampling_factor);
Rx_Buff_dnsamp = [];
for i = 1:upsampling_factor
    Rx_Buff_dnsamp(i,:) = Rx_Buffer(i:upsampling_factor:end);
end

% [Data_freq_off Upchirp_ind] = UC_location_corr(Rx_Buffer(1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind);
[Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer(1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind,16,0.1);
Upchirp_ind
[Data_freq_off, Peak, Upchirp_ind] = dnsamp_buff(Rx_Buff_dnsamp,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre);
[Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre,Fs,Peak);
% [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives_old(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre,Fs,Peak);


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
        if( min(abs(indices(i) - temp(:,1))) > 10 )
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
num_data_sym = 28;%43;%63;%39;
% load('SF8/SF8_symbols.mat');
% load('SF8/SF8_symbols.mat');
% sym = sym(1:num_data_sym);
sym = 0;
for j =1:size(Pream_ind,1)
    j
%     [demod_sym(j,:) sym_peak(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
    [demod_sym(j,:),~] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
    demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j) - 2,N);
end
% 7*mean(Peak_amp(:,3))
%%
load('SF8/SF8_symbols.mat');
num_TX = size(Pream_ind,1);
sym = mod(sym(1:num_data_sym)-2,N);
sym = repmat(sym,num_TX,1);
ser = ((num_data_sym * num_TX) - sum(sum(demod_sym == sym)))/(num_data_sym * num_TX)
pkt_ser = num_data_sym - sum(demod_sym == sym,2)
%%
i = 2;
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
%%
figure
up_thresh = (Peak_amp(:,1) + 4)';
low_thresh = (Peak_amp(:,1) - 4)';
h = boxplot(sym_peak')                        % sample boxplot with grouping
set(h,'LineWidth',2)
hAx=gca;                                   % retrieve the axes handle
xtk=hAx.XTick;                             % and the xtick values to plot() at...
hold on
hL=plot(xtk,up_thresh.*ones(1,size(up_thresh,2)),'k*','MarkerSize',20);     % draw a line on top..
hL=plot(xtk,low_thresh.*ones(1,size(up_thresh,2)),'k*','MarkerSize',20);
set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
xlabel('Packet #','FontSize',30);
ylabel('Peak Height Variation','FontSize',30);
grid minor
ylim([0 30])

%%          UPSAMPLED  DATA

%%      DC correlations
ind = DC_location_correlation(Rx_Buffer,upsampling_factor*N,DC_upsamp(1:upsampling_factor*N),16,0.2);

temp = [];
indices = [zeros(1,floor(num_DC)); ind];

for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) > 3 )
%         temp = [temp; indices(i,:)];
%     end
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 32 )
            temp = [temp; indices(i,:)];
        end
    end
end
DC_ind = temp

%%      UC correlations

[Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer,upsampling_factor*N,num_preamble,num_sync,num_DC,num_data_sym,DC_upsamp,BW,SF,Fs,DC_ind,64,0.1);

temp = [];
indices = [zeros(1,floor(num_preamble)); Upchirp_ind];

for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) > 3 )
%         temp = [temp; indices(i,:)];
%     end
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 32 )
            temp = [temp; indices(i,:)];
        end
    end
end
Upchirp_ind = temp

%%
[Data_freq_off] = dnsamp_buff(Rx_Buffer,Upchirp_ind,num_preamble,num_sync,num_DC,upsampling_factor*N,DC_pre_upsamp);
[Preamble_ind, bin_offsets, Data_out] = filter_false_postives(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,upsampling_factor*N,DC_pre_upsamp,Fs);


temp = [];
temp_data = [];
indices = [zeros(1,num_preamble); Preamble_ind];
Data = [zeros(1,size(Data_out,2)); Data_out];

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
    else
        if( min(abs(indices(i) - temp(:,1))) > 3 )
            temp = [temp; indices(i,:)];
            temp_data = [temp_data; Data(i,:)];
        end
    end
end
Pream_ind = temp
Data_out = temp_data;

%%  Data Demodulations
demod_sym = [];
num_data_sym = 42;%39;
for j =1:size(Pream_ind,1)
    [demod_sym(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind);
%     demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j),N);
end
%%
load('SF8_symbols.mat');
num_TX = size(Pream_ind,1);
sym = repmat(sym,num_TX,1);
((num_data_sym * num_TX) - sum(sum(demod_sym == sym)))/(num_data_sym * num_TX)

%%
i = 10659;%4159;%
buff = Rx_Buffer(1:8:end);
en = sqrt(sum(abs(buff(i:i+N-1)).^2)/length(buff(i:i+N-1)));
FFT = abs(fft(buff(i:i+N-1).*conj(DC)))./sqrt(N);
pnts = find(FFT > en);
plot(FFT,'linewidth',2)
hold on
plot(ones(1,N).*(en),'linewidth',2)
plot(ones(1,N).*(2*en),'linewidth',2)
set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
xlabel('Freq-bin','FontSize',30);
ylabel('Amplidue','FontSize',30);
title('FFT','FontSize',30);

%%      DC correlations
ind = DC_location_correlation(Rx_Buffer(1:upsampling_factor:end),N,DC(1:N),16,0.2);

temp = [];
indices = [zeros(1,floor(num_DC)); ind];

for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) > 3 )
%         temp = [temp; indices(i,:)];
%     end
    if(isempty(temp))
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 3 )
            temp = [temp; indices(i,:)];
        end
    end
end
DC_ind = temp
