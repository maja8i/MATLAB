function YUV = rgbToYuvOffset(RGB)
% RGBTOYUV Convert RGB to YUV
%
% Converts RGB values in [0,255]^3 to YUV values at full scale.
%   Thus Y takes values in the full range [0,255]
%   and U, V may take values outside this range.
%
% RGB = Nx3 input array (uint8, single, double, etc.)
%
% YUV = Nx3 output array (single)

RGB = single(RGB); % convert to floating point if necessary

A = single([    0.299,      0.587,      0.114;...
                -0.14713,   -0.28886,   0.436;...
                0.615,      -0.51499,  -0.10001]);

YUV = RGB * A';
YUV(:,2:3) = YUV(:,2:3) + single(128); % add offset
end

