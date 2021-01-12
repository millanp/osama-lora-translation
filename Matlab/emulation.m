close all;
% clear all;
clc
clearvars -except cursor_info cursor_info1 tx_det

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
% num_samples = pkt_len * N;

DC = conj(sym_to_data_ang([1],N));

UC_SFD = sym_to_data_upsampled([1 1 1],N,Fs,BW);
UC_SFD = UC_SFD(1:2.25*8*N);
DC_upsamp = conj(sym_to_data_upsampled([1],N,Fs,BW));

DC_pre = conj( sym_to_data_ang(preamble_sym  .*  ones(1,num_preamble),N) );
DC_pre_upsamp = conj( sym_to_data_upsampled(preamble_sym  .*  ones(1,num_preamble),N,Fs,BW) );
sync = sym_to_data_ang([9 17],N);
start_frame = [sym_to_data_ang(ones(1,num_preamble),N) sync conj(UC_SFD)];

% fi_1 = fopen('adafruit_SF7');
fi_1 = fopen('SF8/1_tx');
% fi_1 = fopen('pkt_collision');
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

subplot(211)
plot(real(x_1))
subplot(212)
plot(real(x_1_dnsamp))
%%
Packet_number = 1;%input('Which Packet you want to proceed with : ');
% ST_LORA = pkts(:,1);
% ED_LORA = pkts(:,2);
ST_LORA = 1531201;%877456;%cursor_info.DataIndex;
ED_LORA = 1647052;%1380495;%cursor_info1.DataIndex;

template = x_1(ST_LORA(Packet_number):ED_LORA(Packet_number));
num_samples = length(template);
% plot(real(Rx_Buffer))

% plot(real(template(1,1:8:end)))

tx_det = [];
for m = 1:10
    for n = 1:20 
        num_TX = m;
%%
        offset = randi(num_samples,num_TX,1);
        % o = sort(offset)+1 + ((num_preamble + num_sync)*N);
        offset = sort(offset);
        % o(:,2) = o(:,1) + N;
        a = 0.3;
        b = 1;
        temp = zeros(num_TX,num_samples + max(offset(1:num_TX)));
        for i = 1:num_TX
            r = (b-a).*rand(1,1) + a;
        %     temp(i,offset(i)+1 : offset(i)+num_samples) = template;%r*template;
            temp(i,offset(i)+1 : offset(i)+num_samples) = r*template;
        end
        Rx_Buffer_new = sum(temp,1);
        Rx_Buffer = [Rx_Buffer_new zeros(1,3*8*N)];
        Rx_Buffer = awgn(Rx_Buffer,20,'measured');
%         plot(real(Rx_Buffer))
%%      DC correlations
ind = DC_location_correlation(Rx_Buffer(1:upsampling_factor:end),N,DC(1:N),16,0.15);

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
        
        if(size(Pream_ind,1) == 0)
            continue;
        end
        %%  Index check
        true_ind = floor(offset/8) + 165;
        count = 0;
        for i = 1:length(true_ind)
            if(sum(true_ind(i) == Pream_ind(:,1)) || sum(true_ind(i)+1 == Pream_ind(:,1)) || sum(true_ind(i)-1 == Pream_ind(:,1)))
                count = count + 1;
            end
        end
        tx_det(m,n) = count;
    end
end
%%
T = (sum(tx_det,2)/20)';
% T = T ./[1:10] .*100;
b = bar([1:10],T);
ylim([0 10])
set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
% xlabel('x - value','FontSize',30);
xlabel('# of TX collision','FontSize',30);
ylabel('no. of TX located correctly','FontSize',30);
% view(0,90);
set(gca,'YDir','normal');
set(gcf,'Color','w');
% legend('2 TX collision','3 TX collision','4 TX collision')
% legend('SNR = 20dB','SNR = 15dB','SNR = 10dB','SNR = 5dB','SNR = 0dB')
% legend('SNR = 20dB')
% legend('FTrack','CT','CT + offset-filter','CT + PWR-filter','CT + offset-filter + PWR-filter')
% legend('SNR = 20dB','SNR = 15dB','SNR = 10dB','SNR = 5dB')
title('New DC location performance','FontSize',30);
% pbaspect([1 1 1])
grid minor