function [data] = sym_to_data_ang(symbol,N,Fs,BW)
%SYM_TO_DATA Summary of this function goes here
%   Detailed explanation goes here
    data = [];
    accumulator = 0;
    f0 = -BW/2;
    fn = BW/2;
    
    F = f0:BW/(N*(Fs/BW)):fn;
    F = F(1:end-1);
    phase = 0;
    
    for j = symbol
        temp = [];
        c = (j - 1)*(Fs/BW) + 1;
        for i = 1:N*(Fs/BW)
            phase = phase + ((2*pi*F(c))/Fs);
            polar_radius = 1;

            [x, y] = pol2cart(phase, polar_radius);

            temp(i) = complex(x, y);
            c = c+1;
            if(c == N*(Fs/BW))
                c = 1;
            end
        end
        data = [data temp];
    end
end

