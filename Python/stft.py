# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES, WORKING!!!
import numpy as np
import math

def stft(Rx_Buffer,N,DC,upsamp,dis):
#STFT Summary of this function goes here
#   Detailed explanation goes here
#     close all;
    Spec = np.zeros((N,Rx_Buffer.shape[0]))
    # Spec_n = zeros(N,length(Rx_Buffer));
    buff = np.concatenate([Rx_Buffer, np.zeros(N-1)]) # OKAY
    if(upsamp):
        for i in range(Rx_Buffer.shape[1]):
#             Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));#-(i-1)
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC.conj())) / math.sqrt(N),-round( (i)/8 ))#-(i-1)
    else:
        for i in range(len(Rx_Buffer)):
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N), -(i))#
#             Spec(:,i) = abs(fft(buff(i:i+N-1).*DC))./sqrt(N);#-(i-1)
        #     Spec_n(:,i) = circshift(fft(buff(i:i+N-1).*DC.*WinFun)./sqrt(N),-(i-1));
    return Spec

