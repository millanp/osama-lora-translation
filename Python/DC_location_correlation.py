# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES
# NOTE FULLY WORKING!

import numpy as np
import math
from util import length
from get_4_max import get_4_max

def DC_location_correlation(Rx_Buffer,N,DC,pnts_threshold,corr_threshold):
#DC_LOCATION Summary of this function goes here
#   Detailed explanation goes here
#   Detecting downchirp

    Downchirp_ind = []
    Cross_Corr = []

    for i in range(length(Rx_Buffer) - length(DC) - 1):
        Cross_Corr.append(sum(Rx_Buffer[ i : i + (N) ] * DC.conj()) \
                / math.sqrt(sum( Rx_Buffer[ i : i + (N) ] * Rx_Buffer[ i : i + (N) ].conj() ) *
                sum( DC * DC.conj())))
    # figure
    # plot(abs(Cross_Corr))
    # # [~,Downchirp_ind] = max(Cross_Corr);
    # 
    # set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
    # title('Correlation with single Downchirp','FontSize',30);
    # xlabel('Samples','FontSize',30);
    # ylabel('Amp.','FontSize',30);
    # ylim([0 1])
    # keyboard
    # NOTE: OKAY CROSS CORR
    n_samp_array = []
    peak_ind_prev = np.array([])
    Cross_Corr = np.array(Cross_Corr)
    # print('cross corr: ' + str(Cross_Corr))
    for i in range(math.floor(length(Cross_Corr)/N)):
        wind = np.abs(Cross_Corr[i*N : (i+1) * N + 1]) #NOTE OKAY
        # NOTE THRESHOLDS OKAY
        peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold)
        # print('peak ind curr: ' + str(peak_ind_curr[:10]))
    #     if(i == 24)
    #             keyboard
    #         end
        if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):
            
            for j in range(length(peak_ind_curr)):
                for k in range(length(peak_ind_prev)):
                    
                    if(peak_ind_curr[j] == peak_ind_prev[k]):
    #                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                        n_samp_array += [peak_ind_prev[k]+((i-1)*N) + 1, peak_ind_curr[j]+((i)*N) + 1]
        
        peak_ind_prev = peak_ind_curr

    n_samp_array = np.array(n_samp_array) # NOTE OKAY
    for i in range(length(n_samp_array)):
        c = 0
        ind_arr = [n_samp_array[i], n_samp_array[i] + N]
        
        for j in range(len(ind_arr)):
            c = c + np.sum( n_samp_array == ind_arr[j] )

        if( c >= 2 ):
            Downchirp_ind += [ind_arr]
    return np.array(Downchirp_ind)

