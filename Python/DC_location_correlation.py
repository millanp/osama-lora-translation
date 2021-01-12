# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES

import numpy as np
import math
from . import length, get_4_max

def DC_location_correlation(Rx_Buffer,N,DC,pnts_threshold,corr_threshold):
#DC_LOCATION Summary of this function goes here
#   Detailed explanation goes here
#   Detecting downchirp

    Downchirp_ind = []
    Cross_Corr = np.array()

    for i in range(length(Rx_Buffer) - length(DC) - 1):
        Cross_Corr[i] = sum(np.matmul(Rx_Buffer[ i : i + (N) ], DC.conj())) \
                / math.sqrt(sum( np.matmul(Rx_Buffer[ i : i + (N) ], Rx_Buffer[ i : i + (N) ].conj()) ) *
                sum( np.matmul(DC , DC.conj())))
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

    n_samp_array = []
    peak_ind_prev = []
    for i in range(math.floor(length(Cross_Corr)/N)):
        wind = abs(Cross_Corr[i*N + 1 : (i+1) * N + 1])
        peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold)
    #     if(i == 24)
    #             keyboard
    #         end
        if(length(peak_ind_prev) != 0 and length(peak_ind_curr) != 0):
            
            for j in range(length(peak_ind_curr)):
                for k in range(length(peak_ind_prev)):
                    
                    if(peak_ind_curr[j] == peak_ind_prev[k]):
    #                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
    # NOTE: this involves indices... might be off by one
                        n_samp_array = np.concatenate([n_samp_array,  peak_ind_prev[k]+((i)*N), peak_ind_curr[j]+((i + 1)*N)], 1)
        
        peak_ind_prev = peak_ind_curr

    for i in range(length(n_samp_array)):
        c = 0
        ind_arr = [n_samp_array(i), n_samp_array(i) + N]
        
        for j in range(length(ind_arr)):
            c = c + np.sum( n_samp_array == ind_arr[j] )

        if( c >= 2 ):
            Downchirp_ind = np.concatenate([Downchirp_ind, [ind_arr]])
    return Downchirp_ind

