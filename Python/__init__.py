import math
from .sym_to_data_ang import sym_to_data_ang
from .sym_to_data_upsampled import sym_to_data_upsampled
from .get_4_max import get_4_max
from .DC_location_correlation import DC_location_correlation
from .UC_location_corr_DC_based import UC_location_corr_DC_based
from .stft import stft
from .dnsamp_buff import dnsamp_buff
from .Chirplet_Transform import Chirplet_Transform
from .get_bounded_max import get_bounded_max
from .filter_false_postives import filter_false_postives
from .Demod import Demod

def length(arr):
    return max(arr.shape)

def pol2cart(rho, phi):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return(x, y)

import operator as op
from functools import reduce

def nchoosek(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom