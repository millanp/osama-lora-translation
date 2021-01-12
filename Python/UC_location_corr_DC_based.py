# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES (???)
import numpy as np
from . import length, get_4_max
import math

def UC_location_corr_DC_based(Data,N,num_preamble,num_sync,num_DC,num_data_sym,DC,BW,SF,Fs,DC_ind,pnts_threshold,corr_threshold):
    #UC_LOCATION Summary of this function goes here
    #   Detailed explanation goes here

    if(DC_ind.shape[0] == 0):
        return
    pot_pream_ind = np.array([])
    c = 0
    for i in range(DC_ind.shape[0]):
        if(DC_ind[i,0] - ((num_preamble+num_sync)*N) < 1):
            continue
        pot_pream_ind[c,:] = np.arange(DC_ind[i,0] - ((num_preamble + num_sync)*N), DC_ind[i,0]- ((num_sync)*N) + 1, N)
        c = c+1

    Upchirp_ind = np.array([])

    temp_wind = np.array([])
    for j in range(pot_pream_ind.shape[0]):
        if(pot_pream_ind[j,0] - N <= 0):
            continue
        Data_buffer = []
        Data_buffer = Data[pot_pream_ind[j,0] - N : pot_pream_ind[j,-1] + N]
        temp = []
        for i in range(length(Data_buffer) - length(DC)):
            # TODO: why is this i+1???
            temp[i+1] = sum(np.multiply(Data_buffer[i + 1 : i + N], DC[1:N])) \
            / math.sqrt(sum(np.multiply( Data_buffer[i + 1: i + N], Data_buffer[i + 1 : i + N].conj() )) * \
            sum( np.multiply(DC[1:N], DC[1:N].conj())))
        temp_wind[j,:] = temp
    # keyboard
    # figure
    # plot(abs(temp_wind(1,:)));
    # plot(abs(temp_wind(2,:)));

    array_stack = np.array([])
    for m in range(temp_wind.shape[0]):
        
        n_samp_array = []
        peak_ind_prev = []
        for i in range(math.floor(length(temp_wind)/N)):
            
            wind = abs(temp_wind[m,i*N + 1 : (i+1) * N])
            peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold)

            if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):

                for j in range(length(peak_ind_curr)):
                    for k in range(length(peak_ind_prev)):

# TODO TODO TODO TODO: is this an error?
                        if(abs(peak_ind_curr(j)) == abs(peak_ind_prev[k])):
        #                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                            n_samp_array = np.concatenate([n_samp_array,  peak_ind_prev[k]+((i-1)*N)+(pot_pream_ind[m,0]-N-1)], 1)
            peak_ind_prev = peak_ind_curr
        array_stack[m] = n_samp_array

    for m in range(length(array_stack)):
        n_samp_array = np.array(array_stack[m])
        
        for i in range(length(n_samp_array)):
            c = 0
            ind_arr = np.arange(n_samp_array[i] + N, n_samp_array[i] + N + ((num_preamble-2)*N + 1, N))

            for j in range(length(ind_arr)):
                c = c + sum( n_samp_array == ind_arr(j) )
            
            if( c >= 6 ):
                if(length(Upchirp_ind) != 0):
                    if(sum(n_samp_array[i] == Upchirp_ind[:,0]) != 1):
                        Upchirp_ind = np.concatenate([Upchirp_ind, np.concatenate([n_samp_array[i], ind_arr], 1)])
                else:
                    Upchirp_ind = np.concatenate([Upchirp_ind, np.concatenate([n_samp_array[i], ind_arr], 1)])
    temp = np.array([])
    indices = np.concatenate([np.zeros((1,num_preamble)), Upchirp_ind])

    for i in range(1, indices.shape[0]):
        if(length(temp) == 0):
            temp = np.concatenate([temp, indices[i,:]])
        else:
            if( min(abs(indices[i] - temp[:,0])) > 5 ):
                temp = np.concatenate([temp, indices[i,:]])
    Upchirp_ind = temp

    Data_freq_off = 0

    return [Data_freq_off, Upchirp_ind]