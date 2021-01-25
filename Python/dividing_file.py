from sym_to_data_ang import sym_to_data_ang
from sym_to_data_upsampled import sym_to_data_upsampled
from DC_location_correlation import DC_location_correlation
from UC_location_corr_DC_based import UC_location_corr_DC_based
from dnsamp_buff import dnsamp_buff
from filter_false_postives import filter_false_postives
from util import length
from Demod_new import Demod
import numpy as np
import math



num_nodes = 15
num_data_sym = 28
in_file_path = '15tx_0.5lm_1'

for k in [0.5]:
    for iterr in [1]:
        SF = 8
        N = int(2 ** SF)
        BW = 250e3
        Fs = 2e6
        upsampling_factor = int(Fs / BW)
        Ts = 1 / Fs
     
        # SNR = 10;
     
        # Chirplet Transform variables
        sigma = 0.5 #0.94;#0.31623;#0.288675;
        fLevel = N
        WinLen = N
        alpha = - (BW ** 2) / (2 ** SF) #128e6;
     
        # generating LORA pkts
        num_preamble = 8
        num_sync = 2
        preamble_sym = 1
        num_data_sym = 28
        num_DC = 2.25
        pkt_len = num_preamble + num_DC + num_data_sym + 2 # + num_sync
        num_samples = pkt_len * N

        DC = sym_to_data_ang([1], N).conj()
     
        UC_SFD = sym_to_data_upsampled(np.array([1,1,1]), N, Fs, BW)
        UC_SFD = UC_SFD[:int(2.25 * 8 * N)]
        DC_upsamp = sym_to_data_upsampled([1], N, Fs, BW).conj()
     
        DC_pre = sym_to_data_ang(preamble_sym * np.ones(num_preamble), N).conj()
        DC_pre_upsamp = sym_to_data_upsampled(preamble_sym * np.ones((1, num_preamble)), N, Fs, BW).conj()
        sync = sym_to_data_ang(np.array([9, 17]), N)
        start_frame = np.concatenate((sym_to_data_ang(np.ones(num_preamble), N), sync, UC_SFD.conj()))
     
        x_1 = np.fromfile(in_file_path, dtype=np.complex64)
     
        t = np.arange(len(x_1)) / Fs
     
        x_1 = x_1[:int(math.floor(length(x_1) / upsampling_factor) * upsampling_factor)]
     
        x_1_dnsamp = x_1[::int(upsampling_factor)]
     
        chunks = 100
        overlap = 2.5 * upsampling_factor * N
        samp = length(x_1) / chunks
        idx = np.empty((chunks,2))
        for i in range(chunks):
            if (i == 0):
                idx[i, :] = [i * samp, (i+1) * samp]
            else:
                idx[i, :] = [(i * samp) - (overlap), (i+1) * samp]
        file_dur = (length(x_1) / Fs)
        demod_sym_stack = []
        Peaks = []
        for m in range(idx.shape[0]):
            print('Packet ' + str(m+1) + '/' + str(idx.shape[0]))
            ##      DC correlations
            temp_buff = []
            temp_buff = x_1[int(round(idx[m, 0])) : int(idx[m, 1]) + 1]
            ind = DC_location_correlation(temp_buff[::upsampling_factor], N, DC[:N], 16, 0.2)
            if (ind.shape[0] == 0):
                continue
            
            indices = np.concatenate([np.zeros((1, math.floor(num_DC))), ind])
            temp = []
            temp_empty = True
            for i in range(1, indices.shape[0]):
                if (temp_empty):
                    temp.append(indices[i, :])
                    temp_empty = False
                else:
                    if (min(abs(indices[i][0] - np.array(temp)[:, 0])) > 3):
                        temp.append(indices[i, :])
            DC_ind = np.array(temp)
            Rx_Buffer = []
            DC_ind_true = (DC_ind * upsampling_factor) + idx[m, 0]
            for i in range(DC_ind.shape[0]):
                ind1 = int(round(DC_ind_true[i, 0] - (num_preamble + num_sync + 1.5) * upsampling_factor * N))
                ind2 = int(round(DC_ind_true[i, 0] + (num_DC + num_data_sym + 3) * upsampling_factor * N))
                Rx_Buffer.append(x_1[ind1 : ind2 + 1])
            Rx_Buffer = np.array(Rx_Buffer)
            DC_idx = (num_preamble + num_sync + 1.5) * N
            Rx_Buffer = Rx_Buffer[:, :math.floor(Rx_Buffer.shape[1] / upsampling_factor) * upsampling_factor]
         
            ##
            for o in range(Rx_Buffer.shape[0]):
                ##      UC correlations
                Data_freq_off = np.array([])
                Rx_Buff_dnsamp = []
                for i in range(upsampling_factor):
                    Rx_Buff_dnsamp.append(Rx_Buffer[o, i::upsampling_factor])
                Rx_Buff_dnsamp = np.array(Rx_Buff_dnsamp)
             
                [_, Upchirp_ind] = UC_location_corr_DC_based(Rx_Buffer[o, ::upsampling_factor], N, num_preamble, num_sync, num_DC, num_data_sym, DC[:N], BW, SF, Fs, DC_idx, 16, 0.1)
                if (Upchirp_ind.shape[0] == 0):
                    continue
                [Data_freq_off, Peak, Upchirp_ind] = dnsamp_buff(Rx_Buff_dnsamp, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre)
                if (Upchirp_ind.shape[0] == 0):
                    continue
                [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_freq_off, Upchirp_ind, num_preamble, num_sync, num_DC, N, DC_pre, Fs, Peak)
                temp = []
                temp_data = []
                temp_peaks = []
                indices = np.array([np.zeros(num_preamble), Preamble_ind])
                Data = np.array([np.zeros(Data_out.shape[0]), Data_out])
                peaks = np.array([np.zeros(Peak_amp.shape[0]), Peak_amp])
                Peak_amp = []
                for i in range(1, indices.shape[0]):
                    if (len(temp) == 0):
                        temp.append(indices[i, :])
                        temp_data.append(Data[i, :])
                        temp_peaks.append(peaks[i, :])
                    else:
                        if (min(abs(indices[i] - temp[:, 1])) > 10):
                            temp.append(indices[i, :])
                            temp_data.append(Data[i, :])
                            temp_peaks.append(peaks[i, :])
                Pream_ind = np.array(temp)
                Data_out = np.array(temp_data)
                Peak_amp = np.array(temp_peaks)
                ##  Data Demodulations
                demod_sym = []
                sym = 0
                for j in range(Pream_ind.shape[0]):
                    demod_sym.append(Demod(Pream_ind[j, :], Data_out[j, :], BW, SF, Fs, N, num_preamble, num_sync, num_DC, num_data_sym, DC[:N], Pream_ind, Peak_amp[j, :], sym, 7 * np.mean(Peak_amp[:, 2]), m)[0])
                    demod_sym[j] = np.mod(demod_sym[j] + bin_offsets[j] - 2, N)
                demod_sym_stack.append(demod_sym[0])
                Peaks.append(Peak_amp)
        demod_sym_stack = np.array(demod_sym_stack)
        Peaks = np.array(Peaks)
        demod_sym_stack.tofile('symbols.txt')
        Peaks.tofile('peaks.txt')
