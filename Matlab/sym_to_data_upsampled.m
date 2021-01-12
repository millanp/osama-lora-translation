function [data] = sym_to_data_upsampled(symbol,N,Fs,BW)
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
%         phase = -pi + ((j-1)*(2*pi/(N)));
%         phase = 0;
        temp = [];
        c = 1;
        for i = 1:N*(Fs/BW)
            phase = phase + ((2*pi*F(i))/Fs);
            polar_radius = 1;

            [x, y] = pol2cart(phase, polar_radius);

            temp(c) = complex(x, y);
            c = c+1;

        end
        data = [data temp];
    end
end

