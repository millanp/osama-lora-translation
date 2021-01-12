# SYNTAX, FUNCTIONS, MULT, CONCATENATE, INDICES

import numpy as np
import math
from . import pol2cart

def sym_to_data_ang(symbols,N):
#SYM_TO_DATA Summary of this function goes here
#   Detailed explanation goes here
    data = np.array([])
    accumulator = 0
    
    for j in symbols:
        phase = -math.pi + ((j-1)*(2*math.pi/(N)))
#         accumulator = 0;
        temp = []
        for i in range(N):
            accumulator = accumulator + phase
            polar_radius = 1

            [x, y] = pol2cart(accumulator, polar_radius)

            temp[i] = complex(x, y)
#             down_chirp(i) = complex(x, -y);

            phase = phase + (2*math.pi/(N))
        data = np.concatenate([data, temp], 1)
    return data


