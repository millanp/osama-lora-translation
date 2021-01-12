function [Preamble_ind, bin_offsets, Data_out, Peak_amp] = filter_false_postives(Data_stack,Upchirp_ind,num_preamble,num_sync,num_DC,N,DC,Fs,Peak)
%FILTER_FALSE_POSTIVES Summary of this function goes here
%   Detailed explanation goes here

fLevel = N;
WinLen = N;
alpha = 0;

pot_DC_loc = Upchirp_ind(:,1) + ( (num_preamble + num_sync) * N );
Preamble_ind = [];
bin_offsets = [];
Data_out = [];
Peak_amp = [];
for m = 1:size(Upchirp_ind,1)
    data_wind = Data_stack(m,Upchirp_ind(m,1) : Upchirp_ind(m,1) + ((num_preamble )*N) -1);
    c = 1;
    for sigma = [5 ]%0.02]
        close all
        [out,~,~] = Chirplet_Transform(data_wind.*DC,fLevel,WinLen,Fs,alpha,sigma);
%         [temp,~,~] = Chirplet_Transform(data_wind.*DC,fLevel,WinLen,Fs,alpha,sigma);
%         out(c,:,:) = abs(temp);
%         wind(c,:,:) = (abs(temp).^2)./max(max(abs(temp).^2));
%         wind_edge(c,:,:) = edge(abs(temp));
%         spec_plot(temp,N,0)
%         c = c + 1;
    end
%     for i = 1:N+1
%         for j = 1:size(wind,3)
%             absSpec = abs(wind(:,i,j));
%             index = find(min(absSpec) == absSpec);
%             out(i,j) = wind(index(1),i,j);
%         end
%     end
%     spec_plot(out,N,0,0,1);
%     keyboard
    
    temp = [];
    row_ind = [N-4:N+1 1:6];
    count = 1;
    for i = row_ind
        temp(count) = sum(abs(out(i,:)));
        count = count + 1;
    end
    [~,ind] = max(temp);
    pream_peak_ind(m) = row_ind(ind);
    sync1_ind = mod(pream_peak_ind(m) + 8,N);
    sync2_ind = mod(pream_peak_ind(m) + 16,N);
    if(sync1_ind == 0)
        sync1_ind = N;
    end
    if(sync2_ind == 0)
        sync2_ind = N;
    end
    sync_wind = Data_stack(m,Upchirp_ind(m,num_preamble) + N : Upchirp_ind(m,num_preamble) + N + (num_sync*N) - 1);
%     [sync_spec,~,~] = Chirplet_Transform(sync_wind.*DC(1:2*N),fLevel,WinLen,Fs,alpha,0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sync_wind = [zeros(1,N/2) sync_wind zeros(1,N/2)];
%     sync_spec = zeros(N,length(sync_wind) - N);
%     for i = 1:length(sync_wind) - N
% %         en = sqrt(sum(abs(sync_wind(i:i+N-1)).^2)/length(sync_wind(i:i+N-1)));
% %         FFT = circshift(abs(fft(sync_wind(i:i+N-1).*DC(1:N)))./sqrt(N),-(i-1));
%         sync_spec(:,i) = circshift(abs(fft(sync_wind(i:i+N-1).*DC(1:N))),-(i-1-N/2));
% %         pnts = find(FFT > 2.5*en);
% %         if(length(pnts) > 0)
% %             Spec(pnts,i) = FFT(pnts);
% %         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     spec_plot(sync_spec,N,0,0,1);
    
%     sync_threshold = 3*mean(mean(abs(sync_spec)));

    sync_threshold_up = Peak(m,1) + 0.5*Peak(m,1);
    sync_threshold_low = Peak(m,1) - 0.5*Peak(m,1);
    
    sync_word1 = abs(fft(sync_wind(1:N).*DC(1:N)));
    sync_word2 = abs(fft(sync_wind(N+1:end).*DC(1:N)));
%     sync_threshold_up = Peak(m,1) + 7;
%     sync_threshold_low = Peak(m,1) - 7;
    if(sync_threshold_low < (2*sum(sync_word1)/N))
        sync_threshold_low = (2*sum(sync_word1)/N);
    elseif( sync_threshold_low < (2*sum(sync_word2)/N))
        sync_threshold_low = (2*sum(sync_word2)/N);
    end
%     if(sync_threshold_low < 1)
%         sync_threshold_low = 1;
%     end
%     sync_word1 = abs(sync_spec(sync1_ind,1:N));
%     sync_word2 = abs(sync_spec(sync2_ind,N+1:end));

    syn1_pnts = get_bounded_max(sync_word1,sync_threshold_up,sync_threshold_low);%(length(find(sync_threshold < sync_word1))/N);
    syn2_pnts = get_bounded_max(sync_word2,sync_threshold_up,sync_threshold_low);

%     syn1_pnts = (length(find(sync_threshold < sync_word1))/N);
%     syn2_pnts = (length(find(sync_threshold < sync_word2))/N);
    
%     if((length(find(3*mean(mean(out)) < out(pream_peak_ind(m),:)))/size(out,2)) > 0.5 && syn1_pnts > 0.5 && syn2_pnts > 0.5)     %% b == 1 && 
%         Preamble_ind = [Preamble_ind; Upchirp_ind(m,:)];
%         if(pream_peak_ind(m) < N/2)
%             bin_offsets = [bin_offsets 1 + (-mod(pream_peak_ind(m),N))];
%         else
%             bin_offsets = [bin_offsets mod(N+2 - pream_peak_ind(m),N)];
%         end
%         Data_out = [Data_out; Data_stack(m,:)];
%         Peak_amp = [Peak_amp; Peak(m,:)];
% %             count = count+ 1;
%     end

    if(sum(syn1_pnts == sync1_ind) && sum(syn2_pnts == sync2_ind))
        Preamble_ind = [Preamble_ind; Upchirp_ind(m,:)];
        if(pream_peak_ind(m) < N/2)
            bin_offsets = [bin_offsets 1 + (-mod(pream_peak_ind(m),N))];
        else
            bin_offsets = [bin_offsets mod(N+2 - pream_peak_ind(m),N)];
        end
        Data_out = [Data_out; Data_stack(m,:)];
        Peak_amp = [Peak_amp; Peak(m,:)];
%             count = count+ 1;
    end
    
end


end

