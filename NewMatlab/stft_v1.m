function [Spec] = stft_v1(Rx_Buffer,N,DC,upsamp,dis)
%STFT

    % This function produces a spectrogram of LoRa signal using Dechirping
    % operation to get the best frequency Resolution
    Spec = zeros(N,length(Rx_Buffer));
    buff = [Rx_Buffer zeros(1,N-1)];
    if(upsamp)
%         for i = 1:length(Rx_Buffer)
%             Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*conj(DC)))./sqrt(N),-round( (i-1)/8 ));%-(i-1)
%         end
%             spec_plot(abs(Spec),N/8,0,0,0)
    else
        for i = 1:length(Rx_Buffer)
            Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));%\
        end
        if(dis == 1)
%             spec_plot(abs(Spec),N,0,0,0)
        end
    end
end

