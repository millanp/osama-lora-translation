close all;
clear all;
clc

demod_sym_stack = [];
edgeoff = [];

%    clearvars -except iter demod_sym_stack edgeoff
   
   %%
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

load('C:\temp\USRP_data\CT_Matlab\EO\Pkt1');
load('C:\temp\USRP_data\CT_Matlab\EO\Pkt2');
Pkt1_st = 12790;
Pkt2_st = 12789;

delay = 6*N*upsampling_factor;

for m = [delay : delay + (N*upsampling_factor)]
    temp = zeros(2,length(Pkt1) + m);
    temp(1,1:length(Pkt1)) = Pkt1;
    temp(2,m + 1:m + length(Pkt2)) = Pkt2;
    Rx_Buffer = sum(temp,1);
    
% subplot(211)
% plot(real(Pkt1))
% subplot(212)
% plot(real(Rx_Buffer))

    %%  Novel DC Location
    % i = 12159;
    DC_ind = [];
    DC_ind = novel_DC_location(Rx_Buffer(1:8:end),N,DC,num_DC,10);
    DC_ind = sortrows(DC_ind)
    if(size(DC_ind,1) == 0)
        continue;
    end

    %%      UC correlations
    Rx_Buffer = Rx_Buffer(1:floor(length(Rx_Buffer)/upsampling_factor)*upsampling_factor);
    Rx_Buff_dnsamp = [];
    for i = 1:upsampling_factor
        Rx_Buff_dnsamp(i,:) = Rx_Buffer(i:upsampling_factor:end);
    end
    
    Upchirp_ind = [];
    Data_freq_off =[];
    Peak = [];
    bin_offsets =[];
    Data_out =[];
    Peak_amp = [];
    Preamble_ind = [];
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
    if(size(Pream_ind,1)  < 2)
        continue;
    end
    demod_sym = [];
    num_data_sym = 28;%43;%63;%39;
    % load('SF8/SF8_symbols.mat');
    load('EO/sym.mat');
    % sym = sym(1:num_data_sym);
%     sym = 0;
    for j =1:size(Pream_ind,1)
        j
    %     [demod_sym(j,:) sym_peak(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
        [demod_sym(j,:),~] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
        demod_sym(j,:) = mod(demod_sym(j,:) + bin_offsets(j) - 2,N);
    end
    sym = repmat(sym,size(demod_sym,1),1);
    ser = 1 - (sum(sum(demod_sym == sym))/prod(size(sym)));
    
    
edgeoff = [edgeoff; mod(abs(Pream_ind(1,1) - Pream_ind(2,1)),N) ser];
end
% 
% filename = ['edgeoff_6sym/EO_demod'];
% save(filename,'demod_sym_stack');
% filename = ['edgeoff_6sym/EO_offsets'];
% save(filename,'edgeoff');