function [out] = get_4_max(arr,threshold,num_pnts)
%GET_3_MAX Summary of this function goes here
%   Detailed explanation goes here
%     threshold = 50;
    out = [];
    for i = 1:num_pnts
        [a b] = max(arr);
        if(a<threshold)
            return;
        else
            out = [out b];
            arr(b) = 0;
        end
    end
%     [a b] = max(arr);
%     if(a<threshold)
%         return;
%     end
%     out = [out b];
%     arr(b) = 0;
%     [a b] = max(arr);
%     if(a<threshold)
%         return;
%     end
%     out = [out b];
%     arr(b) = 0;
%     [a b] = max(arr);
%     if(a<threshold)
%         return;
%     end
%     out = [out b];
%     arr(b) = 0;
end

