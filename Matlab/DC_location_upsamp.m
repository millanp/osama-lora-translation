function [ind] = DC_location_upsamp(Rx_Buffer,Spec,edge_plot,N,UC,o,UC_SFD)
%DC_LOCATION Summary of this function goes here
%   Detailed explanation goes here
%   Detecting downchirp
% if(~UC)
%     edge_plot = flip(edge_plot,1);
%     Spec = flip(Spec);
% end
edge_range_threshold = 10;
pos_temp = find(abs(edge_plot(N/2,:)) == 1);

pos = [];
for i = 1:length(pos_temp)-1
    if(pos_temp(i + 1) == pos_temp(i) + 1)
        continue;
    else
        if(pos_temp(i + 1) - pos_temp(i) <= edge_range_threshold)
            pos = [pos pos_temp(i) pos_temp(i + 1)];
        end   
    end
end

count = 1;
for j = 1:length(pos)-1
%     if(52391 == pos(j))
%         keyboard
%     end
    if(pos(j + 1) - pos(j) > edge_range_threshold || pos(j + 1) < 1 || pos(j) < 1 || pos(j) > length(Spec) || pos(j+1) > length(Spec))
        continue;
    end
    d = [];
    d_ind = [];
    d = abs(Spec(N/2,pos(j):pos(j+1)));
    d_ind = pos(j):pos(j+1);
    [~,b] = max(d);
    ind_temp(count) = d_ind(b);
    count = count + 1;
end


% filters all points that are N samples apart
temp = [];
for i = ind_temp
    if(sum(i + 8*N == ind_temp) || sum(i + (8*N) + 1 == ind_temp) || sum(i + (8*N) - 1 == ind_temp))
        temp = [temp; i i+(8*N)];
    end
end
keyboard
ind = temp - (8*64);
temp = ind;
for i = 1:size(temp,1)
    data_fft = abs(fft(Rx_Buffer(temp(i,1) : temp(i,1) + (2.25*8*N) - 1) .* UC_SFD));
    plot(data_fft);
end



end

