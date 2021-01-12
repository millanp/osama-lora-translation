function [symbols, sym_peak] = Demod(Pream_ind,Rx_Buffer,BW,SF,Fs,N,num_preamble,num_sync,num_DC,num_data_sym,DC,Pream_frame,Peak_amp,sym,mean_peak_std)
%DEMOD Summary of this function goes here
%   Detailed explanation goes here
% Chirplet Transform variables
% sigma = 0.05;%0.94;%0.94;%0.94;%0.31623;%0.288675;
% 7*(Peak_amp(3))
dis = 0;
edge_range_threshold = 4;
fLevel = N;
WinLen = N;
alpha = 0;%(BW^2)/(2^SF);%128e6; % in Hz/s�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pream_frame = Pream_frame + 1;
Pream_frame(:,num_preamble + 1) = Pream_frame(:,num_preamble) + N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_frame_start = Pream_ind(1) + (num_preamble*N) + (num_DC*N) + (num_sync*N);
Data_frame_end = Data_frame_start + (num_data_sym*N);
% Data_frames = Data_frame_start:N:Data_frame_end;
frame_indices = [((Data_frame_start: N :Data_frame_start+((num_data_sym-1)*N)))' ((Data_frame_start+N-1 : N :Data_frame_start+((num_data_sym)*N)))'];


data_ind = [];
en_arr = [];
sym_peak = [];
for  k = 1:num_data_sym
%     if( sum(k == [ 6 ]))
%         keyboard
%     end
    close all
    %%%%%%%%%%%%%%

    data_wind = Rx_Buffer(frame_indices(k,1):frame_indices(k,2)) .* DC;
    data_wind_next_1 = Rx_Buffer(frame_indices(k,1) + N:frame_indices(k,2) + N) .* DC;
    data_wind_prev_1 = Rx_Buffer(frame_indices(k,1) - N:frame_indices(k,2) - N) .* DC;
    
    data_wind_next_2 = Rx_Buffer(frame_indices(k,1) + 2*N:frame_indices(k,2) + 2*N) .* DC;
    data_wind_prev_2 = Rx_Buffer(frame_indices(k,1) - 2*N:frame_indices(k,2) - 2*N) .* DC;
    temp_d = abs(fft(data_wind,N));
    temp_next_1 = abs(fft(data_wind_next_1,N));
    temp_prev_1 = abs(fft(data_wind_prev_1,N));
    
    temp_next_2 = abs(fft(data_wind_next_2,N));
    temp_prev_2 = abs(fft(data_wind_prev_2,N));
    
    
%     up_thresh = (Peak_amp(1) + 4);
%     low_thresh = (Peak_amp(1) - 4);
    up_thresh = (Peak_amp(1) + 0.25*Peak_amp(1));
    low_thresh = (Peak_amp(1) - 0.25*Peak_amp(1));
%     up_thresh = (Peak_amp(1) + mean_peak_std);
%     low_thresh = (Peak_amp(1) - mean_peak_std);
%     up_thresh = (Peak_amp(1) + 12*(Peak_amp(3)));
%     low_thresh = (Peak_amp(1) - 12*(Peak_amp(3)));
    if(low_thresh < (4*sum(temp_d)/N)) %1
        low_thresh = (4*sum(temp_d)/N);
    end
    
%     plot(temp_d)
%     hold on
%     plot((4*sum(temp_d)/N).*ones(1,N))
%     hold on 
%     plot(up_thresh.*ones(1,N))
%     plot(low_thresh.*ones(1,N))
%     set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
%     title('FFT of Data Window','FontSize',30);
%     xlabel('Freq-bin','FontSize',30);
%     ylabel('Amp.','FontSize',30);
%     keyboard

%     pot_sym = get_4_max(temp_d,max(temp_d)/3,8);
    pot_sym = get_bounded_max(temp_d,up_thresh,low_thresh);
    if(length(pot_sym) == 0)
%         keyboard
    end
%     next_wind_sym_1 = get_4_max(temp_next_1,max(temp_next_1)/3,8);%get_bounded_max(temp_next_1,up_thresh,low_thresh);%
%     next_wind_sym_2 = get_4_max(temp_next_2,max(temp_next_2)/3,8);%get_bounded_max(temp_next_2,up_thresh,low_thresh);%
%     prev_wind_sym_1 = get_4_max(temp_prev_1,max(temp_prev_1)/3,8);%get_bounded_max(temp_prev_1,up_thresh,low_thresh);%
%     prev_wind_sym_2 = get_4_max(temp_prev_2,max(temp_prev_2)/3,8);%get_bounded_max(temp_prev_2,up_thresh,low_thresh);%
%     
    next_wind_sym_1 = get_bounded_max(temp_next_1,up_thresh,low_thresh);%get_4_max(temp_next_1,max(temp_next_1)/3,8);%
    next_wind_sym_2 = get_bounded_max(temp_next_2,up_thresh,low_thresh);%get_4_max(temp_next_2,max(temp_next_2)/3,8);%
    prev_wind_sym_1 = get_bounded_max(temp_prev_1,up_thresh,low_thresh);%get_4_max(temp_prev_1,max(temp_prev_1)/3,8);%
    prev_wind_sym_2 = get_bounded_max(temp_prev_2,up_thresh,low_thresh);%get_4_max(temp_prev_2,max(temp_prev_2)/3,8);%
    
    temp = [];
    for i = 1:length(pot_sym)
        if( (sum(pot_sym(i) == prev_wind_sym_1) && sum(pot_sym(i) == next_wind_sym_1))...
                || (sum(pot_sym(i) == prev_wind_sym_2) && sum(pot_sym(i) == prev_wind_sym_1))...
                || (sum(pot_sym(i) == next_wind_sym_1) && sum(pot_sym(i) == next_wind_sym_2)) )
%             
%         if( (sum(pot_sym(i) == prev_wind_sym_1) || sum(pot_sym(i) == next_wind_sym_1)) )
            
        else
            temp = [temp pot_sym(i)];
        end
    end
    pot_sym = temp;
%     pot_sym = round(pot_sym/2);

    if(length(pot_sym) >= 2)
        r = nchoosek(pot_sym,2);
        freq_diff = abs(r(:,1) - r(:,2));
        freq_diff(find(freq_diff == 1)) = N;
        freq_diff(find(freq_diff == 2)) = N;
        freq_diff;
        min_freq_dif = min(freq_diff);
        sig_f = min_freq_dif / N;
        sig = ((0.05*0.05)/sig_f) + 0.04;
%         if(min(freq_diff) > 20)
%             sig = 0.1;
%         else
%             sig = 0.2;%0.08;
%         end
    else
        sig = 0.1;
    end
    
    % Fractional offset
    
    temp = [];
    for i = 1:length(pot_sym)
        if(sum(pot_sym(i) + 1 == pot_sym) || sum(pot_sym(i) - 1 == pot_sym))
        else
            temp = [temp pot_sym(i)];
        end
    end
    pot_sym = temp;


    %%%%%%%%%%%%%
    if(dis)
        figure;plot(temp);hold on;plot(max(temp)/2.*ones(1,N))
    end
    c = 1;
    wind = [];
    out = [];
    for sigma = [2 sig]
%         close all
        [temp,~,~] = Chirplet_Transform(data_wind,fLevel,WinLen,Fs,alpha,sigma);
         wind(c,:,:) = abs(temp);
%         wind(c,:,:) = (abs(temp).^2)./max(max(abs(temp).^2));
%         wind_edge(c,:,:) = edge(abs(temp));
        if(dis)
            spec_plot(abs(temp),N,0,1,0,0)
        end
%         spec_plot(temp,N,0)
        c = c + 1;
    end
%     spec_plot(reshape(wind(1,:,:),N+1,[]),N,0,1,0,0)
if(dis == 1)
    for i = 1:size(wind,2)
        for j = 1:size(wind,3)
            absSpec = abs(wind(:,i,j));
            index = find(min(absSpec) == absSpec);
            out_temp(i,j) = wind(index(1),i,j);
        end
    end
    spec_plot(out_temp,N,0,1,0,0)
end

    
%     temp = [];
%     for i = 1:length(pot_sym)
%         if(sum(pot_sym(i) == discard_freqs) == 0)
%             temp = [temp pot_sym(i)];
%         end
%     end
%     pot_sym = temp;

    if(length(pot_sym) == 0)
        disp('empty')
        data_ind = [data_ind 0];
%         sym_peak(k) = temp_d(sym(k));
        continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d = [];
    for i = 1:length(pot_sym)
        if(pot_sym(i) > N/2)
            d(i) = pot_sym(i) + 1;
        else
            d(i) = pot_sym(i);
        end
    end
    

    
%     freq_amp = sum(abs(out(d,1:7)),2)/7;
%     freq_amp_end = sum(abs(out(d,end-6:end)),2)/7;

    for i = 1:length(d)
        for j = 1:size(wind,3)
            absSpec = abs(wind(:,d(i),j));
            index = find(min(absSpec) == absSpec);
            out(i,j) = wind(index(1),d(i),j);
        end
    end
    a = [];
    f = [];
    for i = 1:size(out,1)
        a(i,:) = diff(out(i,:));
        f(i) = (length(find(a(i,1:N/2) > 0))*100/N) + (length(find(a(i,N/2+1:end) < 0))*100/N);
    end
%     end
%     freq_amp = abs(out(:,1));
%     freq_amp_end = abs(out(:,end));
% 
%     freq_amp = sum(abs(out(:,1:7)),2)/7;
%     freq_amp_end = sum(abs(out(:,end-6:end)),2)/7;

    freq_amp = sum(abs(out(:,1:14)),2)/14;
    freq_amp_end = sum(abs(out(:,end-13:end)),2)/14;
    
    dif = abs(freq_amp - freq_amp_end);
    [~,b] = min(dif);
%     [~,b] = max(f);
    data_ind = [data_ind pot_sym(b)];
    %sym_peak(k) = temp_d(sym(k));
%     if(k == 18)
%         keyboard
%     end
    
%     [~,s] = max(temp_d);
%     data_ind = [data_ind s];
    
%     en_arr = [en_arr 0 en];
%     s(k) = length(ind);
%     figure;plot([abs(diag(wind,2 - data_ind(2) - 1)).^2; abs(diag(wind,128 + (2 - data_ind(2)) -1)).^2])
    if(dis)
        figure
        plot(abs(out(:,1)))
        hold on
        plot(abs(out(:,end)))
        set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
        % ylim([Freq(1),Freq(end)]);
        view(0,90);
        set(gca,'YDir','normal');
        title('Frequency Spectrum for a particular time instance','FontSize',30);
        xlabel('Freq. / Hz','FontSize',30);
        ylabel('Amplitude','FontSize',30);
        legend('Window Start','Window End')
    end
   
end

symbols = [data_ind; en_arr];

end