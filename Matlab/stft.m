function [Spec] = stft(Rx_Buffer,N,DC,upsamp,dis)
%STFT Summary of this function goes here
%   Detailed explanation goes here
%     close all;
    Spec = zeros(N,length(Rx_Buffer));
    % Spec_n = zeros(N,length(Rx_Buffer));
    buff = [Rx_Buffer zeros(1,N-1)];
    if(upsamp)
        for i = 1:length(Rx_Buffer)
%             Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));%-(i-1)
            Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*conj(DC)))./sqrt(N),-round( (i-1)/8 ));%-(i-1)
        %     Spec_n(:,i) = circshift(fft(buff(i:i+N-1).*DC.*WinFun)./sqrt(N),-(i-1));
        end
        % wind(1,:,:) = abs(Spec);
        % wind(2,:,:) = abs(Spec_n);
        % for i = 1:N
        %     for j = 1:size(wind,3)
        %         out(i,j) = min(abs(wind(:,i,j)));
        % %         absSpec = abs(wind(:,i,j));
        % %         index = find(min(absSpec) == absSpec);
        % %         out(i,j) = wind(index(1),i,j);
        %     end
        % end
        % out = flip(out);
%             Spec = [Spec(N - (N/16)-1 : N,:); Spec(1 : N/16,:)];
            spec_plot(abs(Spec),N/8,0,0,0)
    else
        for i = 1:length(Rx_Buffer)
            Spec(:,i) = circshift(abs(fft(buff(i:i+N-1).*DC))./sqrt(N),-(i-1));%
%             Spec(:,i) = abs(fft(buff(i:i+N-1).*DC))./sqrt(N);%-(i-1)
        %     Spec_n(:,i) = circshift(fft(buff(i:i+N-1).*DC.*WinFun)./sqrt(N),-(i-1));
        end
        if(dis == 1)
            spec_plot(abs(Spec),N,0,0,0)
        end
    end
    
end

