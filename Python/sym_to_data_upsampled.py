# SYNTAX, FUNCTIONS, MULT, CONCATENATE, INDICES(?)

import math
import numpy as np
from util import pol2cart

def sym_to_data_upsampled(symbols,N,Fs,BW):
#SYM_TO_DATA Summary of this function goes here
#   Detailed explanation goes here
    data = []
    accumulator = 0
    f0 = -BW/2
    fn = BW/2
    
    F = np.arange(f0 - 1, fn, BW/(N*(Fs/BW)))#f0:BW/(N*(Fs/BW)):fn
    F = F[:-1]
    phase = 0
    for j in symbols:
#         phase = -pi + ((j-1)*(2*pi/(N)));
#         phase = 0;
        temp = []
        c = 1
        for i in range(N*(int(Fs)//int(BW))): # NOTE switched to int divide
            phase = phase + ((2*math.pi*F[i])/Fs)
            polar_radius = 1

            [x, y] = pol2cart(phase, polar_radius)

            temp.append(complex(x, y))
            c = c+1
        data += temp
    return np.array(data)

