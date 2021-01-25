function [symbols, sym_peak] = Demod(Pream_ind,Rx_Buffer,BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC,Pream_frame,Peak_amp,sym,mean_peak_std,m)
%DEMOD Summary of this function goes here
%   Detailed explanation goes here
% Chirplet Transform variables
% sigma = 0.05;%0.94;%0.94;%0.94;%0.31623;%0.288675;
% 7*(Peak_amp(3))
dis = 0;
wind_min = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pream_frame = Pream_frame;
Pream_frame(:,num_preamble + 1) = Pream_frame(:,num_preamble) + N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(Pream_frame,1)
    frm_st = Pream_frame(i,1) + (num_preamble*N) + (num_DC*N) + (num_sync*N);
    frm_en = frm_st + (num_data_sym*N);
    frm_ind(i,:,:) = [((frm_st: N :frm_st+((num_data_sym-1)*N)))' ((frm_st+N-1 : N :frm_st+((num_data_sym)*N)))'];
end


data_ind = [];
en_arr = [];
sym_peak = [];
    Data_frame_start = Pream_ind(1) + (num_preamble*N) + (num_DC*N) + (num_sync*N);
    Data_frame_end = Data_frame_start + (num_data_sym*N);
    % Data_frames = Data_frame_start:N:Data_frame_end;
    frame_indices = [((Data_frame_start: N :Data_frame_start+((num_data_sym-1)*N)))' ((Data_frame_start+N-1 : N :Data_frame_start+((num_data_sym)*N)))'];
    
    for  k = 1:num_data_sym
        close all;
%         if( sum(k == [ 25 ]))
%             keyboard
%         end
        %%%%%%%%%%%%%%
        %% Find interfering Symbol Boundaries
        ind = [];
        sym_bnd = [];
        for i = 1:size(frm_ind,1)
            if(i == m)
                continue;
            end
%             ind = [ind; reshape(frm_ind(i,intersect(find(frm_ind(i,:,1) > frame_indices(k,1)), find(frm_ind(i,:,1) < frame_indices(k,2))),:),1,[])];
            st = reshape(frm_ind(i,:,1),[],1);
            ed = reshape(frm_ind(i,:,2),[],1);
            sym_bnd = [sym_bnd st(intersect(find(st > frame_indices(k,1)) , find(st < frame_indices(k,2))))];
            sym_bnd = [sym_bnd ed(intersect(find(ed > frame_indices(k,1)) , find(ed < frame_indices(k,2))))];
        end
        
        %% CIC Filtering
        data_wind = Rx_Buffer(frame_indices(k,1):frame_indices(k,2)) .* DC;
        data_fft = abs(fft(data_wind));
        
        sigma = 1;
        WinFun = exp(-(1/(2*(sigma^2)))* linspace(-1,1,N).^2);
        WinFun = WinFun./(sqrt(2*pi)*sigma);
        temp_wind = data_wind .* WinFun;
        
        sym_bnd = mod(sym_bnd - frame_indices(k,1),N);
%         if(sum(sym_bnd < 40) ~= 0)
%             sym_bnd(find(sym_bnd < 40)) = 40;
%             if(sum(sym_bnd > 216) ~= 0)
%                 sym_bnd(find(sym_bnd > 216)) = 216;
%             end
%         end
        intf_wind = [];
        nfft = 4;
        for i = 1:length(sym_bnd)
            buff = zeros(2,nfft*N);
            buff(1,1:sym_bnd(i) - 1) = temp_wind(1:sym_bnd(i) - 1);
            buff(1,:) = abs(fft(buff(1,:),nfft*N))./sqrt(sum(abs(buff(1,:)).^2));
            buff(2,sym_bnd(i):N) = temp_wind(sym_bnd(i):N);
            buff(2,:) = abs(fft(buff(2,:),nfft*N))./sqrt(sum(abs(buff(2,:)).^2));
            intf_wind = [intf_wind; buff];
        end
        intf_wind_min_fft = min(intf_wind,[],1);
        pot_sym_cic = get_4_max(intf_wind_min_fft,4*sum(intf_wind_min_fft)/(nfft*N),nfft*N);
        pot_sym_cic = ceil(pot_sym_cic/nfft);
        %% Power-Filtering
        PwrFctr = 0.5;
        PwrFlr = 4;
        up_thresh = (Peak_amp(1) + PwrFctr*Peak_amp(1));
        low_thresh = (Peak_amp(1) - PwrFctr*Peak_amp(1));
        if(low_thresh < (PwrFlr*sum(data_fft)/N)) %1
            low_thresh = (PwrFlr*sum(data_fft)/N);
        end
        pot_sym_pf = get_bounded_max(data_fft,up_thresh,low_thresh);
        %% Filtering Preamble of interfering Packets
        data_wind_next_1 = Rx_Buffer(frame_indices(k,1) + N:frame_indices(k,2) + N) .* DC;
        data_wind_prev_1 = Rx_Buffer(frame_indices(k,1) - N:frame_indices(k,2) - N) .* DC;
        data_wind_next_2 = Rx_Buffer(frame_indices(k,1) + 2*N:frame_indices(k,2) + 2*N) .* DC;
        data_wind_prev_2 = Rx_Buffer(frame_indices(k,1) - 2*N:frame_indices(k,2) - 2*N) .* DC;
        temp_next_1 = abs(fft(data_wind_next_1,N));
        temp_prev_1 = abs(fft(data_wind_prev_1,N));
        temp_next_2 = abs(fft(data_wind_next_2,N));
        temp_prev_2 = abs(fft(data_wind_prev_2,N));        
        next_wind_sym_1 = get_bounded_max(temp_next_1,up_thresh,low_thresh);%get_4_max(temp_next_1,max(temp_next_1)/3,8);%
        next_wind_sym_2 = get_bounded_max(temp_next_2,up_thresh,low_thresh);%get_4_max(temp_next_2,max(temp_next_2)/3,8);%
        prev_wind_sym_1 = get_bounded_max(temp_prev_1,up_thresh,low_thresh);%get_4_max(temp_prev_1,max(temp_prev_1)/3,8);%
        prev_wind_sym_2 = get_bounded_max(temp_prev_2,up_thresh,low_thresh);%get_4_max(temp_prev_2,max(temp_prev_2)/3,8);%

        temp = [];
        for i = 1:length(pot_sym_pf)
            if( (sum(pot_sym_pf(i) == prev_wind_sym_1) && sum(pot_sym_pf(i) == next_wind_sym_1))...
                    || (sum(pot_sym_pf(i) == prev_wind_sym_2) && sum(pot_sym_pf(i) == prev_wind_sym_1))...
                    || (sum(pot_sym_pf(i) == next_wind_sym_1) && sum(pot_sym_pf(i) == next_wind_sym_2)) )
    %             
    %         if( (sum(pot_sym(i) == prev_wind_sym_1) || sum(pot_sym(i) == next_wind_sym_1)) )

            else
                temp = [temp pot_sym_pf(i)];
            end
        end
        pot_sym_pf = temp;

        %%  Freq. Offset Filtering
%         temp = [];
%         for i = 1:length(pot_sym_pf)
%             if(sum(pot_sym_pf(i) + 1 == pot_sym_pf) || sum(pot_sym_pf(i) - 1 == pot_sym_pf))
%             else
%                 temp = [temp pot_sym_pf(i)];
%             end
%         end
%         pot_sym = temp;
        %%  Choir Module
        npnt = 16;
        data_fft_npnt = abs(fft(data_wind,npnt*N));
        FO_thresh = 0.25;
        sym_FO = [];
        temp = [];
        for i = 1:length(pot_sym_pf)
            ind = [];
            if(pot_sym_pf(i) == 1)
                ind = [(N*npnt) - (npnt/2) + 1 : (N*npnt) (((pot_sym_pf(i)-1) * npnt) + 1) + (npnt/2) : N*npnt];
            else
                ind = (((pot_sym_pf(i)-1) * npnt) + 1) - (npnt/2) : (((pot_sym_pf(i)-1) * npnt) + 1) + (npnt/2);
            end
            [~,a] = max(data_fft_npnt(ind));
            sym_FO = [sym_FO abs(a - ((npnt/2)+1))/npnt];
            if(sym_FO(end) < FO_thresh)
                temp = [temp pot_sym_pf(i)];
            end
        end
%         sym(k)
%         pot_sym_pf - 2
%         sym_FO
        pot_sym = temp;
        
        
        %%
        b = [];
        if(length(sym_bnd) == 0)
            if(length(pot_sym) == 0)
                [~,symbols(k)] = max(data_fft);
            else
                dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
                [~,b] = min(dist);
                symbols(k) = pot_sym(b);
            end
        else
            fin_sym = intersect(pot_sym_cic,pot_sym);
            %%  Final Decision
            if(length(fin_sym) == 0)
                if(length(pot_sym_cic) == 0 && length(pot_sym) ~= 0)
                    dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
                    [~,b] = min(dist);
                    symbols(k) = pot_sym(b);
                elseif(length(pot_sym) == 0 && length(pot_sym_cic) ~= 0)
%                     [~,b] = max(intf_wind_min_fft(nfft.*pot_sym_cic));
%                     symbols(k) = pot_sym_cic(b);
                    sdev = std(intf_wind(:,nfft.*pot_sym_cic),1);
                    [~,b] = min(sdev);
                    symbols(k) = pot_sym_cic(b);
                elseif(length(pot_sym) == 0 && length(pot_sym_cic) == 0)
                    [~,symbols(k)] = max(data_fft);
                else
                    dist = abs(data_fft(pot_sym) - (up_thresh + low_thresh)/2);
                    [~,b] = min(dist);
                    symbols(k) = pot_sym(b);
                end
            else
                %% Max Peak Decision
%                 [~,b] = max(intf_wind_min_fft(fin_sym));
%                 symbols(k) = fin_sym(b);
                %%  Standard Deviation Decision
%                 sdev = std(intf_wind(:,nfft.*fin_sym),1);
%                 [~,b] = min(sdev);
%                 symbols(k) = fin_sym(b);
                %% Dynamic Sigma
%                 avg_pnts = 13;
%                 WinFun1 = exp(-(1/(2*(sig^2)))* linspace(-1,1,N).^2);
%                 WinFun1 = WinFun1./(sqrt(2*pi)*sig);
%                 G_wind1 = data_wind .* WinFun1;
%                 
%                 WinFun2 = exp(-(1/(2*(2^2)))* linspace(-1,1,N).^2);
%                 WinFun2 = WinFun2./(sqrt(2*pi)*2);
%                 G_wind2 = data_wind .* WinFun2;
%                 for i = 0:avg_pnts
%                     Spec_l(1:N/2 + i,i+1) = G_wind1(1:N/2 + i);
%                     Spec_l(N/2 - (avg_pnts - i):N,i+1 + (avg_pnts + 1)) = G_wind1(N/2 - (avg_pnts - i):N);
%                     
%                     Spec_u(1:N/2 + i,i+1) = G_wind2(1:N/2 + i);
%                     Spec_u(N/2 - (avg_pnts - i):N,i+1 + (avg_pnts + 1)) = G_wind2(N/2 - (avg_pnts - i):N);
%                 end
%                 Sp(1,:,:) = fft(Spec_l);
%                 Sp(2,:,:) = fft(Spec_u);
%                 Spec = reshape(min(Sp,[],1),256,[]);
% %                 spec_plot(abs(Spec),N,0,0,0,0);
% 
%                 freq_amp = sum(abs(Spec(fin_sym,1:(avg_pnts + 1))),2)/(avg_pnts + 1);
%                 freq_amp_end = sum(abs(Spec(fin_sym,end-avg_pnts:end)),2)/(avg_pnts + 1);
% 
%                 dif = abs(freq_amp - freq_amp_end);
%                 [~,b] = min(dif);
%                 symbols(k) = fin_sym(b);
                %%  Stft
                avg_pnts = 10;
                G_wind1 = data_wind;% .* WinFun1;
                Spec = [];
                for i = 0:avg_pnts
                    Spec(1:N/2 + i,i+1) = G_wind1(1:N/2 + i);
                    Spec(N/2 - (avg_pnts - i):N,i+1 + (avg_pnts + 1)) = G_wind1(N/2 - (avg_pnts - i):N);
                end
                Spec = fft(Spec);
%                 spec_plot(abs(Spec),N,0,0,0,0);

%                 freq_amp = sum(abs(Spec(fin_sym,1:(avg_pnts+1))),2)/(avg_pnts+1);
%                 freq_amp_end = sum(abs(Spec(fin_sym,end-avg_pnts:end)),2)/(avg_pnts+1);
% 
                freq_amp = min(abs(Spec(fin_sym,1:(avg_pnts+1))),[],2);
                freq_amp_end = min(abs(Spec(fin_sym,end-avg_pnts:end)),[],2);

%                 freq_amp = abs(Spec(fin_sym,1));
%                 freq_amp_end = abs(Spec(fin_sym,end));

                dif = abs(freq_amp - freq_amp_end);
                [~,b] = min(dif);
                symbols(k) = fin_sym(b);
                %%

            end
        end

        end


end
