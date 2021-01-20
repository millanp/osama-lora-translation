# SYNTAX, FUNCTIONS, MULTI, CONCAT, IDX
import numpy as np

def get_4_max(arr,threshold,num_pnts):
#GET_3_MAX Summary of this function goes here
#   Detailed explanation goes here
#     threshold = 50;
    out = []
    for i in range(num_pnts):
        # print('arr: ' + str(arr))
        [a, b] = arr.max(0), arr.argmax(0)
        if(a<threshold):
            return np.array(out) # NOTE: previously returned None
        else:
            out.append(b)
            arr[b] = 0
#     [a b] = max(arr);
#     if(a<threshold)
#         return;
#     end
#     out = [out b];
#     arr(b) = 0;
#     [a b] = max(arr);
#     if(a<threshold)
#         return;
#     end
#     out = [out b];
#     arr(b) = 0;
#     [a b] = max(arr);
#     if(a<threshold)
#         return;
#     end
#     out = [out b];
#     arr(b) = 0;
    return np.array(out)

