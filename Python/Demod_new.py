import numpy as np
from get_bounded_max import get_bounded_max
from Chirplet_Transform import Chirplet_Transform
from util import length, nchoosek
from get_4_max import get_4_max
import math

def Demod(Pream_ind,Rx_Buffer,BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC,Pream_frame,Peak_amp,sym,mean_peak_std,m):
    # Chirplet Transform variables
    dis = 0
    wind_min = 35
    ###########################################################################
    Pream_frame = np.append(Pream_frame, Pream_frame[:,num_preamble-1] + N)
    ###########################################################################
    frm_ind = []
    for i in range(1):
        frm_st = Pream_frame[i] + (num_preamble*N) + (num_DC*N) + (num_sync*N)
        frm_en = frm_st + (num_data_sym*N)
        frm_ind.append([np.arange(frm_st, frm_st+((num_data_sym-1)*N)+1, N), \
                        np.arange(frm_st+N-1, frm_st+((num_data_sym)*N)+1, N)])
    frm_ind = np.array(frm_ind) - 1 # subtract 1 for index conversion
    data_ind = []
    en_arr = []
    sym_peak = []
    Data_frame_start = Pream_ind[0] + (num_preamble*N) + (num_DC*N) + (num_sync*N)
    Data_frame_end = Data_frame_start + (num_data_sym*N)
    frame_indices = np.array([(np.arange(Data_frame_start, Data_frame_start+((num_data_sym-1)*N)+1, N)),
                                    (np.arange(Data_frame_start+N-1, Data_frame_start+((num_data_sym)*N)+1, N))]).T.astype('int') - 1 # subtract 1 for zero indexing
    data_ind = []
    sym_peak = []
    symbols = []
    for  k in range(num_data_sym):
        ## Find interfering Symbol Boundaries
        ind = []
        sym_bnd = []
        for i in range(len(frm_ind)):
            if(i == m):
                continue
            st = frm_ind[i][0]
            ed = frm_ind[i][1]
            newst = st[np.intersect1d((st > frame_indices[k,0]).nonzero() , (st < frame_indices[k,1]).nonzero())]
            newed = ed[np.intersect1d((ed > frame_indices[k,0]).nonzero() , (ed < frame_indices[k,1]).nonzero())]
            if len(newst) != 0:
                sym_bnd.append(newst)
            if len(newed) != 0:
                sym_bnd.append(newed)
        ## CIC Filtering
        data_wind = Rx_Buffer[frame_indices[k,0]:frame_indices[k,1]+1] * DC
        data_fft = abs(np.fft.fft(data_wind))
        sigma = 1
        WinFun = np.exp(-(1/(2*(sigma^2)))* np.linspace(-1,1,N) ** 2)
        WinFun = WinFun / (math.sqrt(2*math.pi)*sigma)
        temp_wind = data_wind * WinFun
        
        sym_bnd = np.mod(sym_bnd - frame_indices[k,0],N)
        intf_wind = []
        nfft = 4
        for i in range(length(sym_bnd)):
            buff = np.zeros((2,nfft*N))
            buff[0,1:sym_bnd(i) - 1] = temp_wind[:sym_bnd[i] - 1]
            buff[0,:] = abs(np.fft.fft(buff[0,:],nfft*N)) / math.sqrt(sum(abs(buff[0,:]) ** 2)) # div OKAY
            buff[1,sym_bnd[i]:N] = temp_wind[sym_bnd[i]:N]
            buff[1,:] = abs(np.fft.fft(buff[1,:],nfft*N)) / math.sqrt(sum(abs(buff[1,:]) ** 2))
            intf_wind.append(buff)
        intf_wind = np.array(intf_wind)
        intf_wind_min_fft = min(intf_wind,[],1)
        pot_sym_cic = get_4_max(intf_wind_min_fft,4*sum(intf_wind_min_fft)/(nfft*N),nfft*N)
        pot_sym_cic = np.ceil(pot_sym_cic/nfft)
        ## Power-Filtering
        PwrFctr = 0.5
        PwrFlr = 4
        up_thresh = (Peak_amp[0] + PwrFctr*Peak_amp[0])
        low_thresh = (Peak_amp[0] - PwrFctr*Peak_amp[0])
        if(low_thresh < (PwrFlr*sum(data_fft)/N)): #1
            low_thresh = (PwrFlr*sum(data_fft)/N)
        pot_sym_pf = get_bounded_max(data_fft,up_thresh,low_thresh)
        ## Filtering Preamble of interfering Packets
        data_wind_next_1 = Rx_Buffer[frame_indices[k,0] + N : frame_indices[k,1] + N+1] * DC
        data_wind_prev_1 = Rx_Buffer[frame_indices[k,0] - N : frame_indices[k,1] - N+1] * DC
        data_wind_next_2 = Rx_Buffer[frame_indices[k,0] + 2*N : frame_indices[k,1] + 2*N+1] * DC
        data_wind_prev_2 = Rx_Buffer[frame_indices[k,0] - 2*N : frame_indices[k,1] - 2*N+1] * DC
        temp_next_1 = abs(np.fft.fft(data_wind_next_1,N))
        temp_prev_1 = abs(np.fft.fft(data_wind_prev_1,N))
        temp_next_2 = abs(np.fft.fft(data_wind_next_2,N))
        temp_prev_2 = abs(np.fft.fft(data_wind_prev_2,N))    
        next_wind_sym_1 = get_bounded_max(temp_next_1,up_thresh,low_thresh)
        next_wind_sym_2 = get_bounded_max(temp_next_2,up_thresh,low_thresh)
        prev_wind_sym_1 = get_bounded_max(temp_prev_1,up_thresh,low_thresh)
        prev_wind_sym_2 = get_bounded_max(temp_prev_2,up_thresh,low_thresh)

        temp = []
        for i in range(length(pot_sym_pf)):
            if( (sum(pot_sym_pf[i] == prev_wind_sym_1) and sum(pot_sym_pf[i] == next_wind_sym_1))\
                    or (sum(pot_sym_pf[i] == prev_wind_sym_2) and sum(pot_sym_pf[i] == prev_wind_sym_1))\
                    or (sum(pot_sym_pf[i] == next_wind_sym_1) and sum(pot_sym_pf[i] == next_wind_sym_2)) ):
                pass                
            else:
                temp.append(pot_sym_pf[i])
        pot_sym_pf = np.array(temp) + 1 # convert from indices to data
        ##  Freq. Offset Filtering
        ##  Choir Module
        npnt = 16
        data_fft_npnt = abs(np.fft.fft(data_wind,npnt*N))
        FO_thresh = 0.25
        sym_FO = []
        temp = []
        for i in range(length(pot_sym_pf)):
            ind = []
            if(pot_sym_pf[i] == 1):
                ind = np.concatenate([np.arange((N*npnt) - (npnt/2) + 1, (N*npnt)), np.arange((((pot_sym_pf[i]-1) * npnt) + 1) + (npnt/2), N*npnt)])
            else:
                ind = np.arange((((pot_sym_pf[i]-1) * npnt) + 1) - (npnt/2), (((pot_sym_pf[i]-1) * npnt) + 1) + (npnt/2) + 1)
            # ind OKAY
            ind = ind.astype('int') - 1 # convert back from data to indices
            a = np.argmax(data_fft_npnt[ind])
            sym_FO.append(abs(a - ((npnt/2)+1))/npnt)
            if(sym_FO[-1] < FO_thresh):
                temp.append(pot_sym_pf[i])
        pot_sym = np.array(temp)
        
        b = []
        if(length(sym_bnd) == 0):
            if(len(pot_sym) == 0):
                symbols.append(np.argmax(data_fft) + 1)
                print('1 1: ', np.argmax(data_fft) + 1)
            else:
                dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2) # subtract 1 from pot_sym to convert to indices

                b = np.argmin(dist)
                symbols.append(pot_sym[b])
                if symbols[-1] == 249:
                    print('1 2: ', pot_sym[b])
        else:
            fin_sym = np.intersect1d(pot_sym_cic,pot_sym) + 1 # convert indices to data
            ##  Final Decision
            if(length(fin_sym) == 0):
                if(length(pot_sym_cic) == 0 and length(pot_sym) != 0):
                    dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2)
                    b = np.argmin(dist)
                    symbols.append(pot_sym[b])
                    
                elif(length(pot_sym) == 0 and length(pot_sym_cic) != 0):
                    sdev = np.std(intf_wind[:,nfft * pot_sym_cic])
                    b = np.argmin(sdev)
                    symbols.append(pot_sym_cic[b])
                    
                elif(length(pot_sym) == 0 and length(pot_sym_cic) == 0):
                    symbols.append(np.argmax(data_fft) + 1) # convert index to data
                    
                else:
                    dist = abs(data_fft[pot_sym - 1] - (up_thresh + low_thresh)/2)
                    b = np.argmin(dist)
                    symbols.append(pot_sym[b])
                    
            else:
                ##  Stft
                avg_pnts = 10
                G_wind1 = data_wind# .* WinFun1;
                Spec = []
                for i in range(avg_pnts+1):
                    Spec[:N/2 + i,i+1] = G_wind1[:N/2 + i]
                    Spec[N/2 - (avg_pnts - i):N,i+1 + (avg_pnts + 1)] = G_wind1[N/2 - (avg_pnts - i):N]
                Spec = np.fft.fft(Spec)
                freq_amp = min(abs(Spec[fin_sym,:(avg_pnts+1)]),[],2)
                freq_amp_end = min(abs(Spec[fin_sym,-avg_pnts:]),[],2)

                dif = abs(freq_amp - freq_amp_end)
                b = np.argmin(dif)
                symbols.append(fin_sym[b])
                

    return [symbols, sym_peak]

