function [Downchirp_ind] = DC_location(Spec,edge_plot,N,UC,o)
%DC_LOCATION Summary of this function goes here
%   Detailed explanation goes here
%   Detecting downchirp
if(~UC)
    edge_plot = flip(edge_plot,1);
    Spec = flip(Spec);
end
edge_range_threshold = 5;
pos = find(abs(edge_plot(N/2,:)) == 1) - (N/2-1);%(N/2 - 1)

temp = [];
for i = 1:length(Spec) - (N)
    temp(:,i) = diag(Spec,i-1);
end
Spec = temp;

count = 1;
for j = 1:length(pos)-1
    if(pos(j + 1) - pos(j) > edge_range_threshold || pos(j + 1) < 1 || pos(j) < 1 || pos(j) > length(Spec) || pos(j+1) > length(Spec))
        continue;
    end
    d = [];
    for i = pos(j):pos(j+1)
        d(i) = sum(abs(Spec(:,i)).^2);
        if(i == pos(j) || i == pos(j+1))
%             find(abs(Spec(:,i)).^2 > sum(abs(Spec(:,i)).^2)/129)
           sum_ones(i) = sum(abs(Spec(:,i)).^2)/129 * 100;
        end
    end
    [~,ind(count)] = max(d);
    count = count+1;
end 

% Dont go inside- just a debug point
if(0)
    count = 0;
    missed = [];
    off = o(:);
    for i = 1:length(off)
        count = count + sum(off(i) == ind);
        if(sum(off(i) == ind) == 0)
            missed = [missed off(i)];
        end
    end
    count
    missed
end

% filters all points that are N samples apart
temp = [];
for i = ind
    if(sum(i + N == ind) || sum(i + N + 1== ind) || sum(i + N - 1 == ind))
        temp = [temp; i i+N];
    end
end
ind = temp;

% temp = zeros(size(Spec));
% for i = ind(1:end-1)
%     vec = abs(diag(Spec,i-1)).';
%     for j = 1:N+1
%         temp(j,i+j-1) = vec(j);
%     end
% end
% Spec = temp;
if(0)
% if(~UC)
    for i = 1:size(ind,1)
        temp = [];
        for  j = ind(i,1)-1:ind(i,1)+1
            temp = [temp; Spec(:,j).'];
        end
        thresh(1,1) = sum(sum(abs(temp).^2))/prod(size(temp));
            
%         s(i) = (max(abs(diag(temp,ind(i)-1)).^2) - min(abs(diag(temp,ind(i)-1)).^2))/2;%mean(abs(diag(temp,ind(i)-1)));
%         threshold = (sum(sum(abs(Spec).^2))/(N*N))/2;%(max(abs(diag(Spec,ind(i)-1)).^2) + min(abs(diag(Spec,ind(i)-1)).^2))/2;mean((abs(diag(Spec,ind(i)-1)).^2).');;%%mean(abs(diag(temp,ind(i)-1)).^2) - std(abs(diag(temp,ind(i)-1)).^2);%

%         threshold = mean((abs(diag(Spec,ind(i)-1)).^2).');
%         threshold1 = std((abs(diag(Spec,ind(i)-1)).^2).');
%         plot(abs(diag(Spec,ind(i)-1)).^2);
%         plot(threshold.*ones(1,N));
%         s(i) = (length(find((abs(diag(Spec,ind(i)-1)).^2) > threshold)) + length(find((abs(diag(Spec,ind(i+1)-1)).^2) > threshold)))/(2*129)*100;%find(abs(diag(Spec,ind(i)-1)).^2 - threshold);

%         s(i) = length(find(abs(Spec(:,ind(i))).^2 > threshold))/(N+1)*100;
        
        a = ((abs(Spec(:,ind(i,:))).^2).');%./max((abs(diag(Spec,ind(i)-1)).^2));
        a(1,:) = smoothdata(a(1,:),'movmean',10);
        a(2,:) = smoothdata(a(2,:),'movmean',10);
%         plot(a)
        a_diff(1,:) = diff(a(1,:));
        a_diff(2,:) = diff(a(2,:));
        pnts(i) = (length(find(a_diff(1,1:N/2) > 0)) + length(find(a_diff(1,N/2+1:end) < 0)) + length(find(a_diff(2,1:N/2) > 0)) + length(find(a_diff(2,N/2+1:end) < 0)))/(2*(N+1)) * 100;
%         e(i) = sum((abs(diag(Spec,ind(i)-1)).^2))/(N+1);
%             a = ((abs(diag(Spec,ind(i)-1)).^2).')./max((abs(diag(Spec,ind(i)-1)).^2));
%             [p,S] = polyfit([1:129],a,2);
%             [y,delta] = polyval(p,[1:129],S);
%             v(i) = sum(abs(a - y));
%             hold on
%             plot(a)
%             plot(y)
    end
    if(length(ind)~=0)
%         figure
    %     stem(ind,s)
    
%         figure
%         stem(ind(:,1),pnts,'linewidth',3,'MarkerSize',10);
%         xlabel('obtained indices','FontSize',30);
%         ylabel('%age of points above and below zero','FontSize',30);
%         set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
        
        
    %     ind = ind(find(s > 20)); %mean(s)
    %     ind = ind(find(v < 12));
        ind = ind(find(pnts > 50),:);
    end
        % last check to get all indices in pairs
        Downchirp_ind = ind + 1;%ind;   % if not taken fftshift
%         for i = min(ind): max(ind)
%             if(sum(ind == i) * sum(ind == i+N))
%                 Downchirp_ind = [Downchirp_ind; i i+N];
%     %             disp('Downchirp found');
%             end
%         end
    end
Downchirp_ind = ind + 1;
end

