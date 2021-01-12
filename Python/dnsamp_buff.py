# SYNTAX, FUNCTIONS, MULT (ALL EW), CONCAT, INDICES(??)
import numpy as np
from . import length, stft
import math

def dnsamp_buff(Data_stack,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC):
#DNSAMP_BUFF Summary of this function goes here
#   Detailed explanation goes here

    Up_ind = []
    peak_amp = []
    Data_buff = []
    n_pnt = 16
    peak_stats = np.array([])
    for k in range(Upchirp_ind.shape[0]):
        if(Upchirp_ind[k,0] - N <= 0):
            continue
        inn = []
        for m in range(Data_stack.shape[0]):
    #         pream_fft = abs(fft(Data_stack(m,Upchirp_ind(k,1) : Upchirp_ind(k,2) - 1) .* DC(1:N)));
    #         [~,pream_bin] = max(pream_fft);
    #         dnchirp_fft = abs(fft(Data_stack(m,Upchirp_ind(k,1) + (num_preamble+num_sync)*N : Upchirp_ind(k,1) + ((num_preamble+num_sync)*N) + N - 1) .* conj(DC(1:N))));
    #         [~,DC_bin] = max(dnchirp_fft);
    #         figure
    #         plot(pream_fft)
    #         figure
    #         plot(dnchirp_fft)
    #         keyboard
                data_wind = np.array([])
                data_fft = np.array([])
                freq_off = np.array([])
                ind_temp = np.concatenate([np.arange(1,5*n_pnt +1), np.arange((N*n_pnt)-(4*n_pnt), (N*n_pnt)+1)], 1)
                c = np.array([]) # NOTE: I added this. where does this variable "c" come from???
                Data_freq_off = np.array([]) # NOTE: I added this. where does this variable come from???
                for j in range(num_preamble):
                    data_wind = Data_stack[m,Upchirp_ind[k,0] : Upchirp_ind[k,0] + (num_preamble*N) -1]
                    data_fft[j,:] = abs(np.fft.fft(data_wind[(j-1)*N + 1:j*N] * DC[1:N],n_pnt*N))
                    
    #                 plot(data_fft(j,:))
                    [_,c[j]] = data_fft[j,ind_temp].max(0)
                    c[j] = ind_temp[c[j]]
                    if(c[j] > (n_pnt*N)/2):
                        freq_off = np.concatenate([freq_off, ( (N*n_pnt) - c[j] ) / n_pnt], 1)
                    else:
                        freq_off = np.concatenate([freq_off, -1*( c(j) - 1 ) / n_pnt], 1)
                        
                freq_off = sum( freq_off[2:7] ) / (num_preamble - 2)
                # NOTE: where did this variable come from???
                Data_freq_off[m,:] = Data_stack[m,:] * math.exp( (1j*2*math.pi*(freq_off / N)) * np.arange(1, length(Data_stack[m,:]) + 1) )
                
                # clear data_wind data_fft ind_temp
                ind_temp = np.concatenate([range(1,6), range(N-4, N+1)], 1)
                a = []
                for j in range(num_preamble):
                    data_wind = Data_freq_off[m,Upchirp_ind[k,0] : Upchirp_ind[k,0] + (num_preamble*N)]
                    data_fft[j,:] = abs(np.fft.fft(data_wind[(j-1)*N + range(1,j*N+1)] * DC[1:N+1],N))
    #                 plot(data_fft(j,:))
                    [a[j],c[j]] = data_fft[j,ind_temp].max(0)
                    c[j] = ind_temp[c[j]]
                peak_stats[k,m,0] = np.mean(a)
    #             if(std(a) < 0.1)
    #                 peak_stats(k,m,2) = 0.1;
    #             elseif(std(a) > 0.3)
    #                 peak_stats(k,m,2) = 0.3;
    #             else
                peak_stats[k,m,1] = np.var(a)#(max(a) - mean(a)) + (mean(a) - min(a));
                peak_stats[k,m,2] = np.std(a)
    #             end
                
            
    #         v = 15;
    #         temp_wind = [];
    #         for i = -v:v
    #             temp_wind(i + v + 1) = sum(Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1) .* DC)...
    #             ./ sqrt(sum( Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1) .* conj(Data_freq_off(m,Upchirp_ind(k,1) + i : Upchirp_ind(k,num_preamble) + N + i - 1)) ) .* ...
    #             sum( DC .* conj(DC)));
    #         end
    #         figure
    #         plot([-v:v],abs(temp_wind))
    #         keyboard

                Spec = stft(Data_freq_off[m,Upchirp_ind[k,0] - N:Upchirp_ind[k,-1] + N - 1 - N],N,DC[1:N],0,0)
                temp = np.array([])
                freq_track_qual = []
                pream_peak_ind = []
                adj_ind = []
                row_ind = np.concatenate([range(N-5,N+1), range(1,6+1)], 1)
                count = 1
                for i in np.nditer(row_ind):
                    temp[count] = np.sum(np.abs(Spec[i,:]))
                    count = count + 1
                [_,ind] = temp.max(0)
                pream_peak_ind = row_ind(ind)
                adj_ind = np.concatenate([np.mod(pream_peak_ind-1,N), np.mod(pream_peak_ind+1,N)], 1)
                if(sum(adj_ind == 0) == 1):
                    adj_ind[(adj_ind == 0).nonzero()] = N
                freq_track_qual = ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind(1),:])) ) + ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind(2),:])) )
    #             keyboard
                inn = np.concatenate([inn, freq_track_qual], 1)


    #         in = [in max(abs(temp_wind))];

        [_,b] = inn.max(0)
    #     Data_buff(k,:) = Data_freq_off(b,:);
        Data_buff = np.concatenate([Data_buff, Data_freq_off[b,:]])
    #     peak_amp(k,:) = reshape(peak_stats(k,b,:),1,[]);
        peak_amp = np.concatenate([peak_amp, np.reshape(peak_stats[k,b,:],(1,-1))]) # TODO: does -1 work the same way as []?
        Up_ind = np.concatenate([Up_ind, Upchirp_ind[k,:]])
        
    return [Data_buff, peak_amp, Up_ind]
