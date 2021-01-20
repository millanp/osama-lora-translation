# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES, WORKING!!!
import numpy as np
import math
from Chirplet_Transform import Chirplet_Transform
from get_bounded_max import get_bounded_max

def filter_false_postives(Data_stack,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC,Fs,Peak):
    #FILTER_FALSE_POSTIVES Summary of this function goes here
    #   Detailed explanation goes here

    fLevel = N
    WinLen = N
    alpha = 0

    # pot_DC_loc = Upchirp_ind[:,0] + ( (num_preamble + num_sync) * N )
    Preamble_ind = []
    bin_offsets = []
    Data_out = []
    Peak_amp = []
    pream_peak_ind = []
    for m in range(Upchirp_ind.shape[0]):
        # DATA STACK OKAY
        data_wind = Data_stack[m,int(Upchirp_ind[m,0]) - 1 : int(Upchirp_ind[m,0] + ((num_preamble )*N)) - 1]
        # DATA WIND OKAY
        [out,_,_] = Chirplet_Transform(data_wind * DC,fLevel,WinLen,Fs,alpha,5)
        # out OKAY
    #     for i = 1:N+1
    #         for j = 1:size(wind,3)
    #             absSpec = abs(wind(:,i,j));
    #             index = find(min(absSpec) == absSpec);
    #             out(i,j) = wind(index(1),i,j);
    #         end
    #     end
    #     spec_plot(out,N,0,0,1);
    #     keyboard
        
        temp = []
        row_ind = np.concatenate([range(N-5, N+1), range(7)])
        count = 0
        for i in np.nditer(row_ind):
            temp.append(np.sum(np.abs(out[i,:])))
            count = count + 1
        ind = np.argmax(temp)
        pream_peak_ind.append(row_ind[ind])
        sync1_ind = np.mod(pream_peak_ind[m] + 8,N)
        sync2_ind = np.mod(pream_peak_ind[m] + 16,N)
        if(sync1_ind == 0):
            sync1_ind = N
        if(sync2_ind == 0):
            sync2_ind = N
        sync_wind = Data_stack[m,int(Upchirp_ind[m,num_preamble-1] + N - 1): int(Upchirp_ind[m,num_preamble-1] + N + (num_sync*N) - 1)]
        # sync wind OKAY

        sync_threshold_up = Peak[m,0] + 0.5*Peak[m,0]
        sync_threshold_low = Peak[m,0] - 0.5*Peak[m,0]
        
        sync_word1 = abs(np.fft.fft(sync_wind[:N] * DC[:N]))
        sync_word2 = abs(np.fft.fft(sync_wind[N:] * DC[:N]))
    #     sync_threshold_up = Peak(m,1) + 7;
    #     sync_threshold_low = Peak(m,1) - 7;
        if(sync_threshold_low < (2*sum(sync_word1)/N)):
            sync_threshold_low = (2*sum(sync_word1)/N)
        elif( sync_threshold_low < (2*sum(sync_word2)/N)):
            sync_threshold_low = (2*sum(sync_word2)/N)
    #     if(sync_threshold_low < 1)
    #         sync_threshold_low = 1;
    #     end
    #     sync_word1 = abs(sync_spec(sync1_ind,1:N));
    #     sync_word2 = abs(sync_spec(sync2_ind,N+1:end));
        syn1_pnts = get_bounded_max(sync_word1,sync_threshold_up,sync_threshold_low)#(length(find(sync_threshold < sync_word1))/N);
        syn2_pnts = get_bounded_max(sync_word2,sync_threshold_up,sync_threshold_low)


        if(sum(syn1_pnts == sync1_ind) and sum(syn2_pnts == sync2_ind)):
            Preamble_ind = np.concatenate([Preamble_ind, Upchirp_ind[m,:]])
            if(pream_peak_ind[m] < N/2):
                bin_offsets.append(1 + (-np.mod(pream_peak_ind[m]+1,N))) # add one since p_p_i consists of indices
            else:
                print('second branch')
                bin_offsets.append(np.mod(N+2 - pream_peak_ind[m],N))
            Data_out = np.concatenate([Data_out, Data_stack[m,:]])
            Peak_amp = np.concatenate([Peak_amp, Peak[m,:]])
    #             count = count+ 1;
    return [Preamble_ind, bin_offsets, Data_out, Peak_amp]

