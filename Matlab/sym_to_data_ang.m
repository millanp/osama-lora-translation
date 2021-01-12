function [data] = sym_to_data_ang(symbol,N)
%SYM_TO_DATA Summary of this function goes here
%   Detailed explanation goes here
    data = [];
    accumulator = 0;
    
    for j = symbol
        phase = -pi + ((j-1)*(2*pi/(N)));
%         accumulator = 0;
        temp = [];
        for i = 1:N
            accumulator = accumulator + phase;
            polar_radius = 1;

            [x, y] = pol2cart(accumulator, polar_radius);

            temp(i) = complex(x, y);
%             down_chirp(i) = complex(x, -y);

            phase = phase + (2*pi/(N));
        end
        data = [data temp];
    end
end

