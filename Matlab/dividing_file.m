num_nodes = 15;
num_data_sym = 28;
path = '/home/millan/wsp/osama/Matlab';
%%
for k = [0.5]
    for iter = 1
        %%
        % chirp variables
        SF = 8;
        N = 2 ^ SF;
        BW = 250e3;
        Fs = 2e6;
        upsampling_factor = Fs / BW;
        Ts = 1 / Fs;
     
        % SNR = 10;
   
        % Chirplet Transform variables
        sigma = 0.5; %0.94;%0.31623;%0.288675;
        fLevel = N;
        WinLen = N;
        alpha = - (BW ^ 2) / (2 ^ SF); %128e6; % in Hz/s�
     
        % generating LORA pkts
        num_preamble = 8;
        num_sync = 2;
        preamble_sym = 1;
        num_data_sym = 28;
        num_DC = 2.25;
        pkt_len = num_preamble + num_DC + num_data_sym + 2; % + num_sync
        num_samples = pkt_len * N;
        % num_TX = 1;%j;%20;
     
        DC = conj(sym_to_data_ang([1], N));
     
        UC_SFD = sym_to_data_upsampled([1 1 1], N, Fs, BW);
        UC_SFD = UC_SFD(1:2.25 * 8 * N);
        DC_upsamp = conj(sym_to_data_upsampled([1], N, Fs, BW));
     
        DC_pre = conj(sym_to_data_ang(preamble_sym .* ones(1, num_preamble), N));
        DC_pre_upsamp = conj(sym_to_data_upsampled(preamble_sym .* ones(1, num_preamble), N, Fs, BW));
        sync = sym_to_data_ang([9 17], N);
        size(sym_to_data_ang(ones(1, num_preamble), N))
        size(sync)
        size(conj(UC_SFD))
        start_frame = [sym_to_data_ang(ones(1, num_preamble), N) sync conj(UC_SFD)];
     
        fi_1 = fopen('15tx_0.5lm_1');
     
        % fi_1 = fopen([path '\' num2str(num_nodes) 'tx_' num2str(k) 'lm_1']);
        x_inter_1 = fread(fi_1, 'float32');
        fclose(fi_1);
        % if data is complex
        x_1 = x_inter_1(1:2:end) + 1i * x_inter_1(2:2:end);
     
        x_1 = x_1.';
        t = [0:length(x_1) - 1] / Fs;
     
        x_1 = x_1(1:floor(length(x_1) / upsampling_factor) * upsampling_factor);
     
        x_1_dnsamp = x_1(1:upsampling_factor:end);
     
        % subplot(211)
        % plot(real(x_1))
        % subplot(212)
        % plot(real(x_1_dnsamp))
        chunks = 100;
        overlap = 2.5 * upsampling_factor * N;
        samp = length(x_1) / chunks;
        for i = 1:chunks
            if (i == 1)
                idx(i, :) = [((i - 1) * samp + 1) i * samp];
            else
                idx(i, :) = [(((i - 1) * samp) - (overlap) + 1) i * samp];
            end
        end
        file_dur = (length(x_1) / Fs)
     
        demod_sym_stack = [];
        Peaks = [];
        %%
        for m = 9:11%size(idx, 1)
            disp(['PROG m ', num2str(m), '/', num2str(size(idx, 1))]);
            %%      DC correlations
            tic
            temp_buff = [];
            temp_buff = x_1(idx(m, 1) : idx(m, 2));
            disp('ENTERING CORRELATION')
            ind = DC_location_correlation(temp_buff(1:upsampling_factor:end), N, DC(1:N), 16, 0.2);
            if (size(ind, 1) == 0)
                continue;
            end
            temp = [];
            indices = [zeros(1, floor(num_DC)); ind];
            size(indices, 1)
            for i = 1:size(indices, 1)
                %     if(abs(indices(i) - indices(i-1)) > 3 )
                %         temp = [temp; indices(i,:)];
                %     end
                if (isempty(temp))
                    temp = [temp; indices(i, :)];
                else
                    unused_indi = indices(i);
                    unused_temp1 = temp(:, 1);
                    if (min(abs(indices(i) - temp(:, 1))) > 3)
                        temp = [temp; indices(i, :)];
                    end
                end
            end
            DC_ind = temp;
         
            Rx_Buffer = [];
            DC_ind_true = (DC_ind * upsampling_factor) + idx(m, 1);
            for i = 1:size(DC_ind, 1)
                unused_inda = DC_ind_true(i, 1) - (num_preamble + num_sync + 1.5) * upsampling_factor * N;
                unused_indb = DC_ind_true(i, 1) + (num_DC + num_data_sym + 3) * upsampling_factor * N;
                Rx_Buffer(i, :) = x_1(DC_ind_true(i, 1) - (num_preamble + num_sync + 1.5) * upsampling_factor * N : DC_ind_true(i, 1) + (num_DC + num_data_sym + 3) * upsampling_factor * N);
            end
            DC_idx = (num_preamble + num_sync + 1.5) * N;
            Rx_Buffer = Rx_Buffer(:, 1:floor(size(Rx_Buffer, 2) / upsampling_factor) * upsampling_factor);
         
            %%
            for o = 1:size(Rx_Buffer, 1)
                disp(['PROG o ', num2str(o), '/', num2str(size(Rx_Buffer, 1))]);
                %%      UC correlations
             
                Data_freq_off = [];
                Rx_Buff_dnsamp = [];
                for i = 1:upsampling_factor
                    Rx_Buff_dnsamp(i, :) = Rx_Buffer(o, i:upsampling_factor:end);
                end
             
                % [Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer(o,1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind,16,0.1);
                [~, Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer(o, 1:upsampling_factor:end), N, num_preamble, num_sync, num_DC, num_data_sym, DC(1:N), BW, SF, Fs, DC_idx, 16, 0.1);
                Upchirp_ind
                if (size(Upchirp_ind, 1) == 0)
                    continue;
                end
                [Data_freq_off, Peak, Upchirp_ind] = dnsamp_buff(Rx_Buff_dnsamp, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre);
                % Upchirp_ind
                if (size(Upchirp_ind, 1) == 0)
                    continue;
                end
                [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_freq_off, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre, Fs, Peak);
                % [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives_old(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre,Fs,Peak);
             
                temp = [];
                temp_data = [];
                temp_peaks = [];
                indices = [zeros(1, num_preamble); Preamble_ind];
                Data = [zeros(1, size(Data_out, 2)); Data_out];
                peaks = [zeros(1, size(Peak_amp, 2)); Peak_amp];
                clear Peak_amp
                for i = 2:size(indices, 1)
                    disp('looping');
                    %     if(abs(indices(i) - indices(i-1)) <= 3 )
                    %
                    %     else
                    %         temp = [temp; indices(i,:)];
                    %         temp_data = [temp_data; Data(i,:)];
                    %     end
                 
                    if (length(temp) == 0)
                        disp('first if');
                        temp = [temp; indices(i, :)];
                        temp_data = [temp_data; Data(i, :)];
                        temp_peaks = [temp_peaks; peaks(i, :)];
                    else
                        disp('second if');
                        if (min(abs(indices(i) - temp(:, 1))) > 10)
                            temp = [temp; indices(i, :)];
                            temp_data = [temp_data; Data(i, :)];
                            temp_peaks = [temp_peaks; peaks(i, :)];
                        end
                    end
                end
                Pream_ind = temp
                Data_out = temp_data;
                Peak_amp = temp_peaks;
                %%  Data Demodulations
                demod_sym = [];
                k
                % load('SF8/SF8_symbols.mat');
                % load('SF8/SF8_symbols.mat');
                % sym = sym(1:num_data_sym);
                sym = 0;
                for j = 1:size(Pream_ind, 1)
                    j
                    %     [demod_sym(j,:) sym_peak(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
                    [demod_sym(j, :), ~] = Demod(Pream_ind(j, :), Data_out(j, :), BW, SF, Fs, N, num_preamble, num_sync, num_DC, num_data_sym, DC(1:N), Pream_ind, Peak_amp(j, :), sym, 7 * mean(Peak_amp(:, 3)), m);
                    demod_sym(j, :) = mod(demod_sym(j, :) + bin_offsets(j) - 2, N);
                end
                demod_sym_stack = [demod_sym_stack; demod_sym];
                Peaks = [Peaks; Peak_amp];
            end
         
        end
        disp('writing');
        % filename1 = ['C:\Osama\Matlab_lt\25%Peak\' num2str(k) 'tx_lamda' num2str(iter)];
        % filename2 = ['C:\Osama\Matlab_lt\25%Peak\' num2str(k) 'tx_lamda' num2str(iter) '_Peaks'];
        % filename1 = [path '\lamda_' num2str(floor(k * num_nodes)) '_' num2str(iter)];
        % filename2 = [path '\lamda_' num2str(floor(k * num_nodes)) '_' num2str(iter) '_Peaks'];
        save('thing1.txt', 'demod_sym_stack');
        save('thing2.txt', 'Peaks');
    end
 
end