"""
Translation process
1. Correct the syntax of this file
2. Line by line, update the function names etc. from MATLAB to Numpy
- When I encounter a call to another MATLAB function to be translated, move
  into that file and proceed from step 1
3. Fix all matrix multiplications
4. Fix all matrix concatenations
5. Fix all indices (shaky, will be finalized when I run the files)

Later:
- Run the Python file and compare output with the MATLAB program.
- Fix index issues
- Fix array expansion issues
"""

# SYNTAX, 

from . import sym_to_data_ang, sym_to_data_upsampled, length, DC_location_correlation, UC_location_corr_DC_based, dnsamp_buff
import numpy as np
import math



num_nodes = 15
num_data_sym = 28
path = 'C:\Osama\Matlab_lt\indoor_15nodes_30dBspd_out'
for k in [0.5]:
    for iterr in [1]:
        SF = 8
        N = 2 ^ SF
        BW = 250e3
        Fs = 2e6
        upsampling_factor = Fs / BW
        Ts = 1 / Fs
     
        # SNR = 10;
     
        # Chirplet Transform variables
        sigma = 0.5 #0.94;#0.31623;#0.288675;
        fLevel = N
        WinLen = N
        alpha = - (BW ^ 2) / (2 ^ SF) #128e6; # in Hz/s�
     
        # generating LORA pkts
        num_preamble = 8
        num_sync = 2
        preamble_sym = 1
        num_data_sym = 28
        num_DC = 2.25
        pkt_len = num_preamble + num_DC + num_data_sym + 2 # + num_sync
        num_samples = pkt_len * N
        # num_TX = 1;#j;#20;

        DC = sym_to_data_ang([1], N).conj()
     
        UC_SFD = sym_to_data_upsampled(np.concatenate([1,1,1]), N, Fs, BW)
        UC_SFD = UC_SFD[:2.25 * 8 * N]
        DC_upsamp = sym_to_data_upsampled([1], N, Fs, BW).conj()
     
        DC_pre = sym_to_data_ang(preamble_sym * np.ones(1, num_preamble), N).conj()
        DC_pre_upsamp = sym_to_data_upsampled(preamble_sym * np.ones(1, num_preamble), N, Fs, BW).conj()
        sync = sym_to_data_ang(np.concatenate([9, 17]), N)
        start_frame = np.concatenate([sym_to_data_ang(np.ones(1, num_preamble), N), sync, UC_SFD.conj()])
     
        x_1 = np.fromfile('15tx_0.5lm_1', dtype=np.complex64)
     
        # # fi_1 = fopen([path '\' num2str(num_nodes) 'tx_' num2str(k) 'lm_1']);
        # x_inter_1 = fread(fi_1, 'float32');
        # fclose(fi_1);
        # # if data is complex
        # x_1 = x_inter_1(1:2:end) + 1i * x_inter_1(2:2:end);
     
        # x_1 = x_1.';
        t = np.arange(len(x_1)) / Fs
     
        x_1 = x_1[:math.floor(length(x_1) / upsampling_factor) * upsampling_factor]
     
        x_1_dnsamp = x_1[::upsampling_factor]
     
        # subplot(211)
        # plot(real(x_1))
        # subplot(212)
        # plot(real(x_1_dnsamp))
        chunks = 100
        overlap = 2.5 * upsampling_factor * N
        samp = length(x_1) / chunks
        idx = np.array() # TODO NEW
        for i in range(chunks):
            if (i == 1):
                idx[i, :] = np.concatenate([((i - 1) * samp + 1), i * samp])
            else:
                idx[i, :] = np.concatenate([(((i - 1) * samp) - (overlap) + 1), i * samp])
        file_dur = (length(x_1) / Fs)
     
        demod_sym_stack = []
        Peaks = []
        ##
        for m in range(idx.shape[0]):
            ##      DC correlations
            # tic
            temp_buff = []
            temp_buff = x_1[idx[m, 1] : idx[m, 2]]
            ind = DC_location_correlation(temp_buff[::upsampling_factor], N, DC[:N], 16, 0.2)
         
            if (ind.shape[0] == 0):
                continue
            temp = np.array([])
            indices = np.concatenate([np.zeros((1, math.floor(num_DC))), ind])
         
            for i in range(1, indices.shape[0]):
                #     if(abs(indices(i) - indices(i-1)) > 3 )
                #         temp = [temp; indices(i,:)];
                #     end
                if (temp.size == 0):
                    temp = np.concatenate([temp, indices[i, :]])
                else:
                    if (min(abs(indices[i] - temp[:, 1])) > 3):
                        temp = np.concatenate([temp, indices[i, :]])
            DC_ind = temp
         
            Rx_Buffer = np.array([])
            DC_ind_true = (DC_ind * upsampling_factor) + idx(m, 1)
            for i in range(DC_ind.shape[0]):
                Rx_Buffer[i, :] = x_1[DC_ind_true(i, 1) - (num_preamble + num_sync + 1.5) * upsampling_factor * N : DC_ind_true(i, 1) + (num_DC + num_data_sym + 3) * upsampling_factor * N]
            DC_idx = (num_preamble + num_sync + 1.5) * N
            Rx_Buffer = Rx_Buffer[:, :math.floor(Rx_Buffer.shape[1] / upsampling_factor) * upsampling_factor]
         
            ##
            for o in range(Rx_Buffer.shape[0]):
                ##      UC correlations
             
                Data_freq_off = np.array([])
                Rx_Buff_dnsamp = np.array([])
                for i in range(upsampling_factor):
                    Rx_Buff_dnsamp[i, :] = Rx_Buffer[o, i::upsampling_factor]
             
                # [Data_freq_off Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer(o,1:upsampling_factor:end),N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),BW,SF,Fs,DC_ind,16,0.1);
                [_, Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer[o, ::upsampling_factor], N, num_preamble, num_sync, num_DC, num_data_sym, DC[:N], BW, SF, Fs, DC_idx, 16, 0.1)
                print(Upchirp_ind)
                if (Upchirp_ind.shape[0] == 0):
                    continue
                [Data_freq_off, Peak, Upchirp_ind] = dnsamp_buff(Rx_Buff_dnsamp, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre)
                # Upchirp_ind
                if (Upchirp_ind.shape[0] == 0):
                    continue
                # NOTE: drilling in
                [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_freq_off, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre, Fs, Peak)
                # [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives_old(Data_freq_off,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC_pre,Fs,Peak);
             
                temp = []
                temp_data = []
                temp_peaks = []
                indices = np.concatenate([np.zeros((1, num_preamble)), Preamble_ind])
                Data = np.concatenate([np.zeros((1, size(Data_out, 2))), Data_out])
                peaks = np.concatenate([np.zeros((1, size(Peak_amp, 2))), Peak_amp])
                Peak_amp = np.array()# clear Peak_amp
                for i in range(1, indices.shape[0]):
                    #     if(abs(indices(i) - indices(i-1)) <= 3 )
                    #
                    #     else
                    #         temp = [temp; indices(i,:)];
                    #         temp_data = [temp_data; Data(i,:)];
                    #     end
                 
                    if (length(temp) == 0):
                        temp = np.concatenate([temp, indices[i, :]])
                        temp_data = np.concatenate([temp_data, Data[i, :]])
                        temp_peaks = np.concatenate([temp_peaks, peaks[i, :]])
                    else:
                        if (min(abs(indices[i] - temp[:, 1])) > 10):
                            temp = np.concatenate([temp, indices[i, :]])
                            temp_data = np.concatenate([temp_data, Data[i, :]])
                            temp_peaks = np.concatenate([temp_peaks, peaks[i, :]])
                Pream_ind = temp
                Data_out = temp_data
                Peak_amp = temp_peaks
                ##  Data Demodulations
                demod_sym = []
                print(k)
                # load('SF8/SF8_symbols.mat');
                # load('SF8/SF8_symbols.mat');
                # sym = sym(1:num_data_sym);
                sym = 0
                for j in range(size(Pream_ind, 1)):
                    print(j)
                    #     [demod_sym(j,:) sym_peak(j,:)] = Demod(Pream_ind(j,:),Data_out(j,:),BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC(1:N),Pream_ind,Peak_amp(j,:),sym,7*mean(Peak_amp(:,3)));
                    [demod_sym[j, :], _] = Demod(Pream_ind[j, :], Data_out[j, :], BW, SF, Fs, N, num_preamble, num_sync, num_DC, num_data_sym, DC[1:N], Pream_ind, Peak_amp[j, :], sym, 7 * mean(Peak_amp[:, 3]))
                    demod_sym[j, :] = mod(demod_sym[j, :] + bin_offsets(j) - 2, N)
                demod_sym_stack = np.concatenate([demod_sym_stack, demod_sym])
                Peaks = np.concatenate([Peaks, Peak_amp])
     
        # filename1 = ['C:\Osama\Matlab_lt\25#Peak\' num2str(k) 'tx_lamda' num2str(iter)];
        # filename2 = ['C:\Osama\Matlab_lt\25#Peak\' num2str(k) 'tx_lamda' num2str(iter) '_Peaks'];
        # filename1 = [path '\lamda_' num2str(floor(k * num_nodes)) '_' num2str(iter)];
        # filename2 = [path '\lamda_' num2str(floor(k * num_nodes)) '_' num2str(iter) '_Peaks'];
        # save(filename1, 'demod_sym_stack');
        # save(filename2, 'Peaks');
        np.savetxt('out.symbols', demod_sym_stack, dtype=np.complex64)
