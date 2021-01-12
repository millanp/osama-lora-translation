function [Upchirp_ind bin_offsets Data_out] = UC_location(Data,N,UC,o,num_preamble,num_sync,DC,BW,SF,sigma,Fs,pream_peak)
%UC_LOCATION Summary of this function goes here
%   Detailed explanation goes here
fLevel = N;
WinLen = N;
alpha = 0;
% sigma = 0.1;

Upchirp_ind = [];
edge_range_threshold = 4;

c = 1;
if(size(o,1) == 0)
    return;
end
for i = 1:size(o,1)
    if(o(i,1) - ((num_preamble+num_sync)*N) - 1 < 1)
        continue;
    end
    pot_pream_ind(c,:) = o(i,1) - ((num_preamble + num_sync)*N) - 1 : N : o(i,1)- ((num_sync)*N);
    c = c+1;
end

% for i = 1:length(pream_peak)
%     if(pream_peak(i) > N/2)
%         pream_peak_ind(i) = pream_peak(i) + 1;
%     end
% end

% count = 0;
max_fft = [];
max_fft_ind = [];
bin_offsets = [];
Data_out = [];
for k = 1:size(pot_pream_ind,1)
%     for i = -3:3%-N/4 : N/4
        data_wind = Data(k,pot_pream_ind(k,1) + 1:pot_pream_ind(k,num_preamble) + N); %pot_pream_ind(k,num_preamble) + N);
        data_fft = abs(fft(data_wind.*DC));
        c = 1;
        for sigma = [5 0.02]
            close all
            [temp,~,~] = Chirplet_Transform(data_wind.*DC,fLevel,WinLen,Fs,alpha,sigma);
             wind(c,:,:) = abs(temp);
    %         wind(c,:,:) = (abs(temp).^2)./max(max(abs(temp).^2));
    %         wind_edge(c,:,:) = edge(abs(temp));
    %         spec_plot(temp,N,0)
            c = c + 1;
        end

        for i = 1:N+1
            for j = 1:size(wind,3)
                absSpec = abs(wind(:,i,j));
                index = find(min(absSpec) == absSpec);
                out(i,j) = wind(index(1),i,j);
            end
        end
        spec_plot(out,N,0,0,1); %pot_pream_ind(k,num_preamble) + N);
%     keyboard
%     for j = 1:size(wind,3)
%         absSpec = abs(wind(:,1,j));
%         index = find(min(absSpec) == absSpec);
%         out(i,j) = wind(index(1),1,j);
%     end
        
%         
%         max_fft = [max_fft max(temp(i + 3 + 1,:))];
%         max_fft_ind = [max_fft_ind pot_pream_ind(k,1) + 1 + i];
        [~,b] = max(data_fft);
        temp = [];
        row_ind = [N-13:N+1 1:15];
        count = 1;
        for i = row_ind
            temp(count) = sum(abs(out(i,:)));
            count = count + 1;
        end
        [~,ind] = max(temp);
        pream_peak_ind(k) = row_ind(ind);
        sync1_ind = mod(pream_peak_ind(k) + 7,N);
        sync2_ind = mod(pream_peak_ind(k) + 15,N);
        if(sync1_ind == 0)
            sync1_ind = N;
        end
        if(sync2_ind == 0)
            sync2_ind = N;
        end
        sync_wind = Data(k,pot_pream_ind(k,num_preamble) + N + 1 : pot_pream_ind(k,num_preamble) + N + (num_sync*N));
        [sync_spec,~,~] = Chirplet_Transform(sync_wind.*DC(1:2*N),fLevel,WinLen,Fs,alpha,0.5);
        sync_threshold = 3*mean(mean(abs(sync_spec)));
        sync_word1 = abs(sync_spec(sync1_ind,1:N));
        sync_word2 = abs(sync_spec(sync2_ind,N+1:end));
        
        syn1_pnts = (length(find(sync_threshold < sync_word1))/N);
        syn2_pnts = (length(find(sync_threshold < sync_word2))/N);
%         figure;
%         plot(abs(sync_word1))
%         hold on
%         plot(sync_threshold.*ones(1,N))
        
        if((length(find(3*mean(mean(out)) < out(pream_peak_ind(k),:)))/size(out,2)) > 0.5 && syn1_pnts > 0.5 && syn2_pnts > 0.5)     %% b == 1 && 
            Upchirp_ind = [Upchirp_ind; pot_pream_ind(k,:)];
            if(pream_peak_ind(k) < N/2)
                bin_offsets = [bin_offsets 1 + (-mod(pream_peak_ind(k),N))];
            else
                bin_offsets = [bin_offsets mod(N+2 - pream_peak_ind(k),N)];
            end
            Data_out = [Data_out; Data(k,:)];
%             count = count+ 1;
        end
%     end
%         figure
%         plot(abs(out(1,:)))
%         hold on
%         plot(2*mean(mean(out)).*ones(1,size(out,2)))
end
% Upchirp_ind = pot_pream_ind;

% pnts = zeros(size(Upchirp_ind,1),1);
% for i = 1:size(Upchirp_ind,1)
%     a = [];
%     for j = 1:num_preamble
%         temp = abs(diag(Spec,Upchirp_ind(i,j)-1)).^2;
%         UC_avg_energy(i,j) = sum(temp);
%         a(j,:) = diff(smoothdata(temp,'movmean',10));
%         pnts(i) = pnts(i) + length(find(a(j,1:N/2) > 0)) + length(find(a(j,N/2+1 : end) < 0)) ;
%     end
%     
% end
% pnts = pnts./(num_preamble*(N+1)) .* 100;
% % figure
% % stem(Upchirp_ind(:,1),pnts);
% ind = find(pnts > 50);
% Upchirp_ind = Upchirp_ind(ind,:);
end

