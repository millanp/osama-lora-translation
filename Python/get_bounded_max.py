# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES

def get_bounded_max(arr,up_thresh,low_thresh):
    pnts_up = (arr < up_thresh).nonzero()
    pnts_low = (arr > low_thresh).nonzero()
    pnts = pnts_up.intersection(pnts_low)
    return pnts
