# SYNTAX, FUNCTIONS, MULT, 

import numpy as np
from get_bounded_max import get_bounded_max
from Chirplet_Transform import Chirplet_Transform
from util import length, nchoosek
import math

def Demod(Pream_ind,Rx_Buffer,BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC,Pream_frame,Peak_amp,sym,mean_peak_std):
    #DEMOD Summary of this function goes here
    #   Detailed explanation goes here
    # Chirplet Transform variables
    # sigma = 0.05;#0.94;#0.94;#0.94;#0.31623;#0.288675;
    # 7*(Peak_amp(3))
    dis = 0
    edge_range_threshold = 4
    fLevel = N
    WinLen = N
    alpha = 0;#(BW^2)/(2^SF);#128e6; # in Hz/s�
    ###########################################################################
    Pream_frame = Pream_frame + 1
    Pream_frame = np.append(Pream_frame, Pream_frame[:,num_preamble-1] + N)
    ###########################################################################

    Data_frame_start = Pream_ind[0] + (num_preamble*N) + (num_DC*N) + (num_sync*N)
    Data_frame_end = Data_frame_start + (num_data_sym*N)
    # Data_frames = Data_frame_start:N:Data_frame_end;
    frame_indices = np.array([(np.arange(Data_frame_start, Data_frame_start+((num_data_sym-1)*N)+1, N)),
                                    (np.arange(Data_frame_start+N-1, Data_frame_start+((num_data_sym)*N)+1, N))]).T.astype('int') - 1 # subtract 1 for zero indexing
    # frame indices OKAY

    data_ind = []
    sym_peak = []
    for  k in range(num_data_sym):
    #     if( sum(k == [ 6 ]))
    #         keyboard
    #     end
        ##############

        data_wind = Rx_Buffer[frame_indices[k,0]:frame_indices[k,1]+1] * DC
        data_wind_next_1 = Rx_Buffer[frame_indices[k,0] + N : frame_indices[k,1] + N+1] * DC
        data_wind_prev_1 = Rx_Buffer[frame_indices[k,0] - N : frame_indices[k,1] - N+1] * DC
        
        data_wind_next_2 = Rx_Buffer[frame_indices[k,0] + 2*N : frame_indices[k,1] + 2*N+1] * DC
        data_wind_prev_2 = Rx_Buffer[frame_indices[k,0] - 2*N : frame_indices[k,1] - 2*N+1] * DC
        temp_d = abs(np.fft.fft(data_wind,N))
        temp_next_1 = abs(np.fft.fft(data_wind_next_1,N))
        temp_prev_1 = abs(np.fft.fft(data_wind_prev_1,N))
        
        temp_next_2 = abs(np.fft.fft(data_wind_next_2,N))
        temp_prev_2 = abs(np.fft.fft(data_wind_prev_2,N))
        
        
        up_thresh = (Peak_amp[0] + 0.25*Peak_amp[0])
        low_thresh = (Peak_amp[0] - 0.25*Peak_amp[0])
        if(low_thresh < (4*sum(temp_d)/N)): #1
            low_thresh = (4*sum(temp_d)/N)
        
        pot_sym = get_bounded_max(temp_d,up_thresh,low_thresh)
        next_wind_sym_1 = get_bounded_max(temp_next_1,up_thresh,low_thresh)#get_4_max(temp_next_1,max(temp_next_1)/3,8);#
        next_wind_sym_2 = get_bounded_max(temp_next_2,up_thresh,low_thresh)#get_4_max(temp_next_2,max(temp_next_2)/3,8);#
        prev_wind_sym_1 = get_bounded_max(temp_prev_1,up_thresh,low_thresh)#get_4_max(temp_prev_1,max(temp_prev_1)/3,8);#
        prev_wind_sym_2 = get_bounded_max(temp_prev_2,up_thresh,low_thresh)#get_4_max(temp_prev_2,max(temp_prev_2)/3,8);#
        
        temp = []
        for i in range(length(pot_sym)):
            if( (sum(pot_sym[i] == prev_wind_sym_1) and sum(pot_sym[i] == next_wind_sym_1))\
                    or (sum(pot_sym[i] == prev_wind_sym_2) and sum(pot_sym[i] == prev_wind_sym_1))\
                    or (sum(pot_sym[i] == next_wind_sym_1) and sum(pot_sym[i] == next_wind_sym_2)) ):
                pass
    #             
    #         if( (sum(pot_sym(i) == prev_wind_sym_1) || sum(pot_sym(i) == next_wind_sym_1)) )
                
            else:
                temp.append(pot_sym[i])
        pot_sym = np.array(temp)
    #     pot_sym = round(pot_sym/2);

        if(length(pot_sym) >= 2):
            r = nchoosek(pot_sym,2)
            freq_diff = abs(r[:,0] - r[:,1])
            freq_diff[np.nonzero(freq_diff == 1)] = N
            freq_diff[np.nonzero(freq_diff == 2)] = N
            min_freq_dif = min(freq_diff)
            sig_f = min_freq_dif / N
            sig = ((0.05*0.05)/sig_f) + 0.04
        else:
            sig = 0.1
        
        # Fractional offset
        
        temp = []
        for i in range(length(pot_sym)):
            if(sum(pot_sym[i] + 1 == pot_sym) or sum(pot_sym[i] - 1 == pot_sym)):
                pass
            else:
                temp.append(pot_sym[i])
        pot_sym = np.array(temp)


        #############
        c = 0
        wind = []
        
        for sigma in [2, sig]:
    #         close all
            [temp,_,_] = Chirplet_Transform(data_wind,fLevel,WinLen,Fs,alpha,sigma)
            wind.append(abs(temp))
            c = c + 1
        wind = np.array(wind)

        if(length(pot_sym) == 0):
            data_ind.append(0)
    #         sym_peak(k) = temp_d(sym(k));
            continue
        #####################################
        d = []
        for i in range(length(pot_sym)):
            if(pot_sym[i] > N/2):
                d.append(pot_sym[i] + 1)
            else:
                d.append(pot_sym[i])
        
        out = np.empty((len(d), wind.shape[2]))
        
    #     freq_amp = sum(abs(out(d,1:7)),2)/7;
    #     freq_amp_end = sum(abs(out(d,end-6:end)),2)/7;

        for i in range(len(d)):
            for j in range(wind.shape[2]):
                absSpec = abs(wind[:,d[i],j])
                index = np.nonzero(absSpec == min(absSpec))
                out[i,j] = wind[index[0],d[i],j]
        # a = []
        # f = []
        # for i in range(out.shape[0]):
        #     a.append(np.diff(out[i,:], n=1, axis=0))
        #     f.append((length(np.nonzero(a[i,1:N/2] > 0))*100/N) + (length(np.nonzero(a[i,N/2+1:] < 0))*100/N))

        freq_amp = np.sum(abs(out[:,:14]),1)/14
        freq_amp_end = np.sum(abs(out[:,-14:]),1)/14 # changed -13 to -14 due to python indexing
        
        dif = abs(freq_amp - freq_amp_end)
        b = np.argmin(dif)
        # dif OKAY
    #     [~,b] = max(f);
        data_ind.append(pot_sym[b])
    symbols = np.array(data_ind)
    return [symbols, sym_peak]
