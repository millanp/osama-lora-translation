import numpy as np
import math
from util import length
from get_4_max import get_4_max
import time

def DC_location_correlation(Rx_Buffer,N,DC,pnts_threshold,corr_threshold):
#   Detecting downchirp

    Downchirp_ind = []
    Cross_Corr = []

    st = time.perf_counter()
    ###################################################
    # BEGIN SLOW SECTION
    print('num: ', length(Rx_Buffer) - length(DC) - 1)
    # OLD:
    for i in range(length(Rx_Buffer) - length(DC) - 1):
        Cross_Corr.append(np.sum(Rx_Buffer[ i : i + N ] * DC.conj()) \
                / math.sqrt(np.sum( Rx_Buffer[ i : i + N ] * Rx_Buffer[ i : i + (N) ].conj() ) *
                np.sum( DC * DC.conj())))
    # END SLOW SECTION
    ###################################################
    n_samp_array = []
    peak_ind_prev = np.array([])
    Cross_Corr = np.array(Cross_Corr)
    for i in range(math.floor(length(Cross_Corr)/N)):
        wind = np.abs(Cross_Corr[i*N : (i+1) * N + 1])
        peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold)
        if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):
            
            for j in range(length(peak_ind_curr)):
                for k in range(length(peak_ind_prev)):
                    
                    if(peak_ind_curr[j] == peak_ind_prev[k]):
                        n_samp_array += [peak_ind_prev[k]+((i-1)*N) + 1, peak_ind_curr[j]+((i)*N) + 1]
        
        peak_ind_prev = peak_ind_curr

    n_samp_array = np.array(n_samp_array)
    for i in range(length(n_samp_array)):
        c = 0
        ind_arr = [n_samp_array[i], n_samp_array[i] + N]
        
        for j in range(len(ind_arr)):
            c = c + np.sum( n_samp_array == ind_arr[j] )

        if( c >= 2 ):
            Downchirp_ind += [ind_arr]
    return np.array(Downchirp_ind)

