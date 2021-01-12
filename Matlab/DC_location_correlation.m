function [Downchirp_ind] = DC_location_correlation(Rx_Buffer,N,DC,pnts_threshold,corr_threshold)
%DC_LOCATION Summary of this function goes here
%   Detailed explanation goes here
%   Detecting downchirp

Downchirp_ind = [];

for i = 1:length(Rx_Buffer) - length(DC) - 1
    Cross_Corr(i) = sum(Rx_Buffer( i : i + (N) - 1) .* conj(DC))...
            / sqrt(sum( Rx_Buffer( i : i + (N) - 1) .* conj(Rx_Buffer( i : i + (N) - 1)) ) * ...
            sum( DC .* conj(DC)));
end
% figure
% plot(abs(Cross_Corr))
% % [~,Downchirp_ind] = max(Cross_Corr);
% 
% set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
% title('Correlation with single Downchirp','FontSize',30);
% xlabel('Samples','FontSize',30);
% ylabel('Amp.','FontSize',30);
% ylim([0 1])
% keyboard

n_samp_array = [];
peak_ind_prev = [];
for i = 0:floor(length(Cross_Corr)/N)-1
    wind = abs(Cross_Corr(i*N + 1 : (i+1) * N));
    peak_ind_curr = get_4_max(wind,corr_threshold,pnts_threshold);
%     if(i == 24)
%             keyboard
%         end
    if(length(peak_ind_prev) ~= 0 && length(peak_ind_curr) ~= 0)
        
        for j = 1:length(peak_ind_curr)
            for k = 1:length(peak_ind_prev)
                
                if(peak_ind_curr(j) == peak_ind_prev(k))
%                     n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+(i*N)];
                    n_samp_array = [n_samp_array  peak_ind_prev(k)+((i-1)*N) peak_ind_curr(j)+((i)*N)];
                end
                
            end
        end
        
    end
    
    peak_ind_prev = peak_ind_curr;
end

for i = 1:length(n_samp_array)
    c = 0;
    ind_arr = n_samp_array(i) : N : n_samp_array(i) + (N);
    
    for j = 1:length(ind_arr)
        c = c + sum( n_samp_array == ind_arr(j) );
    end

    if( c >= 2 )
        Downchirp_ind = [Downchirp_ind; [ind_arr]];
    end
end

end

