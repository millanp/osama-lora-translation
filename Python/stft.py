import numpy as np
import math

def stft(Rx_Buffer,N,DC,upsamp,dis):
    Spec = np.zeros((N,Rx_Buffer.shape[0]))
    buff = np.concatenate([Rx_Buffer, np.zeros(N-1)])
    if(upsamp):
        for i in range(Rx_Buffer.shape[1]):
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC.conj())) / math.sqrt(N),-round( (i)/8 ))
    else:
        for i in range(len(Rx_Buffer)):
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N), -(i))
    return Spec

