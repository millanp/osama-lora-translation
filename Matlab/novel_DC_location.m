function [DC_ind] = novel_DC_location(buff,N,DC,num_DC,thresh)

Spec = zeros(N,length(buff) - N);
for i = 1:length(buff) - N
    en = sqrt(sum(abs(buff(i:i+N-1)).^2)/length(buff(i:i+N-1)));
    FFT = abs(fft(buff(i:i+N-1).*conj(DC)))./sqrt(N);
    pnts = find(FFT > 2.5*en);
    if(length(pnts) > 0)
        Spec(pnts,i) = FFT(pnts);
    end
end
% spec_plot(abs(Spec),N,0,0,0)

% 

for i = 1:N
    Spec_new(i,:) = circshift(Spec(i,:),-(N - i));
end
% spec_plot(abs(Spec_new),N,0,0,0)
ener = sum(Spec_new,1);
% figure
% plot(ener)

% plot(ener,'linewidth',2)
% hold on
% set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
% xlabel('samples','FontSize',30);
% ylabel('Amplidue','FontSize',30);
% title('New DC Detection based on dynamic threshold','FontSize',30);
tracks = get_4_max(ener,thresh,length(ener));%find(ener > 50);
ind = [];
for i = 1:length(tracks)
    if(sum((tracks(i) + N) == tracks))
        ind = [ind; tracks(i) (tracks(i) + N)];
    end
end

temp = [];
indices = [zeros(1,floor(num_DC)); ind];

for i = 2:size(indices,1)
%     if(abs(indices(i) - indices(i-1)) > 3 )
%         temp = [temp; indices(i,:)];
%     end
    if(length(temp) == 0)
        temp = [temp; indices(i,:)];
    else
        if( min(abs(indices(i) - temp(:,1))) > 10 )
            temp = [temp; indices(i,:)];
        end
    end
end
DC_ind = temp;


end