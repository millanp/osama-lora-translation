# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES, WORKING!!!
import numpy as np
from util import length
from get_4_max import get_4_max
import math

def UC_location_corr_DC_based(Data,N,num_preamble,num_sync,num_DC,num_data_sym,DC,BW,SF,Fs,DC_ind,pnts_threshold,corr_threshold):
    #UC_LOCATION Summary of this function goes here
    #   Detailed explanation goes here

    # if(DC_ind.shape[0] == 0):
    #     return # TODO what does this mean???
    pot_pream_ind = []
    for i in range(1):
        if(DC_ind - ((num_preamble+num_sync)*N) < 1):
            continue
        pot_pream_ind.append(np.arange(DC_ind - ((num_preamble + num_sync)*N), DC_ind- ((num_sync)*N) + 1, N))
    pot_pream_ind = np.array(pot_pream_ind)

    Upchirp_ind = []

    temp_wind = []
    for j in range(pot_pream_ind.shape[0]):
        if(pot_pream_ind[j,0] - N <= 0):
            continue
        Data_buffer = []
        Data_buffer = Data[int(pot_pream_ind[j,0] - N) : int(pot_pream_ind[j,-1] + N)]
        temp = [0+0j]
        for i in range(length(Data_buffer) - length(DC)):
            # TODO: why is this i+1???
            temp.append(sum(np.multiply(Data_buffer[i + 1 : i + N + 1], DC[:N])) \
                / math.sqrt(sum(Data_buffer[i + 1: i + N + 1] * Data_buffer[i + 1 : i + N + 1].conj() ) * \
                sum( DC[:N] * DC[:N].conj())))
        temp_wind.append(temp)
    temp_wind = np.array(temp_wind)
    # keyboard
    # figure
    # plot(abs(temp_wind(1,:)));
    # plot(abs(temp_wind(2,:)));

    array_stack = []
    for m in range(temp_wind.shape[0]):
        
        n_samp_array = []
        peak_ind_prev = np.array([])
        for i in range(math.floor(length(temp_wind)/N)):
            
            wind = abs(temp_wind[m,i*N + 1 : (i+1) * N])
            peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold)

            if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):

                for j in range(length(peak_ind_curr)):
                    for k in range(length(peak_ind_prev)):
                        if(abs(peak_ind_curr[j]) == abs(peak_ind_prev[k])):
        #                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                            n_samp_array.append(peak_ind_prev[k]+((i-1)*N)+(pot_pream_ind[m,0]-N-1) + 3) # TODO why add three? something to do with indices
            peak_ind_prev = peak_ind_curr
        array_stack.append(n_samp_array)
    array_stack = np.array(array_stack) # NOTE OKAY

    for m in range(len(array_stack)):
        n_samp_array = np.array(array_stack[m])
        
        for i in range(length(n_samp_array)):
            c = 0
            ind_arr = np.arange(n_samp_array[i] + N, n_samp_array[i] + N + (num_preamble-2)*N + 1, N)
            for j in range(len(ind_arr)):
                c = c + sum( n_samp_array == ind_arr[j] )
            if( c >= 6 ):
                if(len(Upchirp_ind) != 0):
                    if(sum(np.array(Upchirp_ind)[:,0] == n_samp_array[i]) != 1):
                        Upchirp_ind.append(np.concatenate([[n_samp_array[i]],ind_arr]))
                else:
                    Upchirp_ind.append(np.concatenate([[n_samp_array[i]], ind_arr]))
    Upchirp_ind = np.array(Upchirp_ind)
    temp = []
    indices = np.concatenate([np.zeros((1,num_preamble)), Upchirp_ind])
    # NOTE indices good
    for i in range(1, indices.shape[0]):
        if(len(temp) == 0):
            temp.append(indices[i,:])
        else:
            if( min(abs(indices.flatten(1)[i] - np.array(temp)[:,0])) > 5 ):
                temp.append(indices[i,:])
    temp = np.array(temp)
    Upchirp_ind = temp

    Data_freq_off = 0

    return [Data_freq_off, Upchirp_ind]