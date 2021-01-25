import numpy as np
from util import length
from stft import stft
import math, cmath

def dnsamp_buff(Data_stack,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC):
    Up_ind = []
    peak_amp = []
    Data_buff = []
    n_pnt = 16
    peak_stats = []
    for k in range(Upchirp_ind.shape[0]):
        if(Upchirp_ind[k,0] - N <= 0):
            peak_stats.append([])
            continue
        inn = []
        k_peak_stats = []
        Data_freq_off = []
        for m in range(Data_stack.shape[0]):
            data_wind = []
            data_fft = []
            freq_off = []
            ind_temp = np.concatenate([np.arange(5*n_pnt), np.arange((N*n_pnt)-(4*n_pnt)-1, (N*n_pnt))])
            c = []
            for j in range(num_preamble):
                data_wind = Data_stack[m,int(Upchirp_ind[k,0]) - 1 : int(Upchirp_ind[k,0] + (num_preamble*N) -1)]
                data_fft.append(abs(np.fft.fft(data_wind[((j)*N):((j+1)*N)] * DC[:N],n_pnt*N)))
                
                c.append(data_fft[j][ind_temp].argmax(0))
                c[j] = ind_temp[c[j]] + 1
                if(c[j] > (n_pnt*N)/2):
                    freq_off.append(( (N*n_pnt) - c[j] ) / n_pnt)
                else:
                    freq_off.append(-1*( c[j] - 1 ) / n_pnt)
            freq_off = sum( freq_off[1:7] ) / (num_preamble - 2)
            Data_freq_off.append(Data_stack[m,:] * np.exp( (1j*2*math.pi*(freq_off / N)) * np.arange(1, length(Data_stack[m,:]) + 1) ))
            ind_temp = np.concatenate([range(5), range(N-4, N)])
            a = []
            c = []
            data_wind = []
            data_fft = []
            for j in range(num_preamble):
                data_wind = Data_freq_off[m][int(Upchirp_ind[k,0]) - 1 : int(Upchirp_ind[k,0] + (num_preamble*N)) - 1]
                data_fft.append(abs(np.fft.fft(data_wind[(j)*N : (j+1)*N] * DC[:N],N)))
                [aj,cj] = data_fft[j][ind_temp].max(0), data_fft[j][ind_temp].argmax(0)
                a.append(aj); c.append(cj)

                c[j] = ind_temp[c[j]]
            k_peak_stats.append([np.mean(a), np.var(a, ddof=1), np.std(a, ddof=1)])
            Spec = stft(Data_freq_off[m][int(Upchirp_ind[k,0] - N)-1:int(Upchirp_ind[k,-1] + N - 1 - N)],N,DC[:N],0,0)
            temp = []
            freq_track_qual = []
            pream_peak_ind = []
            adj_ind = []
            row_ind = np.concatenate([range(N-6,N), range(0,6)])
            count = 1
            for i in np.nditer(row_ind):
                temp.append(np.sum(np.abs(Spec[i,:])))
                count = count + 1
            temp = np.array(temp)
            ind = temp.argmax(0)
            pream_peak_ind = row_ind[ind]
            adj_ind = np.array([np.mod(pream_peak_ind-1,N), np.mod(pream_peak_ind+1,N)])
            if(sum(adj_ind == 0) == 1):
                adj_ind[(adj_ind == 0).nonzero()] = N
            freq_track_qual = ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[0],:])) ) + ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[1],:])) )
            inn.append(freq_track_qual)
        inn = np.array(inn)
        peak_stats.append(k_peak_stats)
        Data_freq_off = np.array(Data_freq_off)
        b = inn.argmax(0)
        Data_buff.append(Data_freq_off[b,:])
        peak_amp.append(peak_stats[k][b])
        Up_ind.append(Upchirp_ind[k,:])
    Data_buff = np.array(Data_buff)
    Up_ind = np.array(Up_ind)
    peak_amp = np.array(peak_amp)
    return [Data_buff, peak_amp, Up_ind]
