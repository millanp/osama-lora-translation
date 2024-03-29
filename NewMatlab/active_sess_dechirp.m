function [uplink_wind] = active_sess_dechirp(x_1)
%ACTIVE_SESS_DECHIRP Summary of this function goes here
%   Detailed explanation goes here

SF = param_configs(1);
BW = param_configs(2);
Fs = param_configs(3);
N = 2^SF;
upsampling_factor = Fs/BW;

DC = conj(sym_to_data_ang([1],N));
DC_fft = fft(DC);
DC_upsamp =(ifft([DC_fft(1:N/2) zeros(1,(upsampling_factor-1)*N) DC_fft(N/2 + 1:N)]));

peak_gain = [];
uplink_wind = [];
n = [];
p = [];
last_wind = 0;
win_jump_factor = 3;
front_buf = 6*win_jump_factor;
back_buf = 3*win_jump_factor;
win_jump = floor(N*upsampling_factor/win_jump_factor);
mov_thresh_wind = 1000*win_jump_factor;
mov_thresh = 0;
mov_thresh_rec  = [];
for i = 1:floor(length(x_1)/win_jump) - win_jump_factor%length(DC))
    f = (i-1)*win_jump+1;
    g = (i-1)*win_jump+(N*upsampling_factor);
    wind = x_1((i-1)*win_jump+1 : (i-1)*win_jump+(N*upsampling_factor));
    wind_fft = abs(fft(wind.*DC_upsamp));
    wind_fft = wind_fft([1:N/2 (N/2 + (upsampling_factor-1)*N)+1:(upsampling_factor)*N]);    
    noise_floor = mean(wind_fft);
    n = [n noise_floor];
    fft_peak = max(wind_fft);
    p = [p fft_peak];
    peak_gain = [peak_gain 10*log10(fft_peak/noise_floor)];
    if(i > mov_thresh_wind)
        mov_thresh = 1.3*mean(peak_gain(end - mov_thresh_wind + 1: end));
        if(mov_thresh > 6)
            mov_thresh = 6;
        end
    else
        mov_thresh = 1.3*mean(peak_gain);
        if(mov_thresh > 6)
            mov_thresh = 6;
        end
    end
    mov_thresh_rec = [mov_thresh_rec mov_thresh];
    
    if(peak_gain(end) >= mov_thresh)
        if(i > last_wind)
            if(i-back_buf < 1)
                uplink_wind = [uplink_wind; 1 i+front_buf];
            else
                uplink_wind = [uplink_wind; i-back_buf i+front_buf];
            end
            last_wind = uplink_wind(end,2);
        elseif(i <= last_wind)
            uplink_wind(end,2) = i + front_buf;
            last_wind = uplink_wind(end,2);
        end
    else
        
    end
%     plot(wind_fft)
%     pause(0.1);
end
uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= front_buf + back_buf),:);
uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= (front_buf + back_buf + 1)),:);
uplink_wind = uplink_wind(find(uplink_wind(:,2)-uplink_wind(:,1) ~= (front_buf + back_buf - 1)),:);
temp_link = uplink_wind;
% uplink_wind = uplink_wind.*(win_jump*upsampling_factor);
uplink_wind = uplink_wind.*(win_jump);

% idx = find(uplink_wind(:,2)-uplink_wind(:,1) < ((pkt_len)*N*upsampling_factor));
if(uplink_wind(end,2) > length(x_1))
    uplink_wind(end,2) = length(x_1);
end
% uplink_wind(idx,2) = uplink_wind(idx,1) + ((pkt_len + 2)*N*upsampling_factor);
% peak_gain = movmean(peak_gain,20);

plot(peak_gain,'linewidth',1.5)
hold on
% plot(ones(1,length(peak_gain)).*6)
plot(mov_thresh_rec);
title('Peak Gain Plot');

set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
xlabel('Samples','FontSize',30);
set(gca,'YDir','normal');
set(gcf,'Color','w');
     ylabel('Peak Gain Magnitude','FontSize',30);
grid minor
ylim([0 max(peak_gain)])

% figure
% plot(n)
% hold on
% plot(p)
end

