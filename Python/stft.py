# SYNTAX, FUNCTIONS, MULT, CONCAT, INDICES (?)
import numpy as np
import math

def stft(Rx_Buffer,N,DC,upsamp,dis):
#STFT Summary of this function goes here
#   Detailed explanation goes here
#     close all;
    Spec = np.zeros((N,Rx_Buffer.shape[1]))
    # Spec_n = zeros(N,length(Rx_Buffer));
    buff = np.concatenate([Rx_Buffer, np.zeros((1,N-1))], 1)
    if(upsamp):
        for i in range(Rx_Buffer.shape[1]):
#             Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));#-(i-1)
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC.conj())) / math.sqrt(N),-round( (i)/8 ))#-(i-1)
        #     Spec_n(:,i) = circshift(fft(buff(i:i+N-1).*DC.*WinFun)./sqrt(N),-(i-1));
        # wind(1,:,:) = abs(Spec);
        # wind(2,:,:) = abs(Spec_n);
        # for i = 1:N
        #     for j = 1:size(wind,3)
        #         out(i,j) = min(abs(wind(:,i,j)));
        # #         absSpec = abs(wind(:,i,j));
        # #         index = find(min(absSpec) == absSpec);
        # #         out(i,j) = wind(index(1),i,j);
        #     end
        # end
        # out = flip(out);
#             Spec = [Spec(N - (N/16)-1 : N,:); Spec(1 : N/16,:)];
        # spec_plot(abs(Spec),N/8,0,0,0) # NOTE: removed plot
    else:
        for i in range(Rx_Buffer.shape[1]):
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N), -(i))#
#             Spec(:,i) = abs(fft(buff(i:i+N-1).*DC))./sqrt(N);#-(i-1)
        #     Spec_n(:,i) = circshift(fft(buff(i:i+N-1).*DC.*WinFun)./sqrt(N),-(i-1));
    return Spec

