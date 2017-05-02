function RGB = yuvToRgbOffset(YUV)
% YUVTORGB Convert YUV to RGB
%
% Converts YUV values at full scale to RGB values in [0,255]^3.
%   Thus Y takes values in the full range [0,255]
%   and U, V may take values outside this range.
%
% YUV = Nx3 input array (single)
%
% RGB = Nx3 output array (single)

YUV = single(YUV); % convert to floating point if necessary

A = single([    1,          0,          1.13983;...
                1,          -0.39465,   -0.58060;...
                1,          2.03211,    0]);
            
YUV(:,2:3) = YUV(:,2:3) - single(128); % remove offset

RGB = YUV * A';

end

