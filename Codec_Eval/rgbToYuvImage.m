function YUV = rgbToYuvImage(RGB)
% Convert RGB to YUV
%
% Converts RGB values in [0,255]^3 to YUV values at full scale.
%   Thus Y takes values in the full range [0,255]
%   and U, V may take values outside this range.
%
% RGB = MxNx3 input array (uint8, single, double, etc.), where MxN are the
%       dimensions of the input image
%
% YUV = MxNx3 output array (single), where YUV(:, :, 1) corresponds to Y
%       values; YUV(:, :, 2) corresponds to U values; and YUV(:, :, 3)
%       corresponds to V values

RGB = single(RGB); % convert to floating point if necessary

A = single([    0.299,      0.587,      0.114;...
                -0.14713,   -0.28886,   0.436;...
                0.615,      -0.51499,  -0.10001]);

YUV(:, :, 1) = RGB(:, :, 1).* A(1, 1) + RGB(:, :, 2).* A(1, 2) + RGB(:, :, 3).* A(1, 3);
YUV(:, :, 2) = RGB(:, :, 1).* A(2, 1) + RGB(:, :, 2).* A(2, 2) + RGB(:, :, 3).* A(2, 3);
YUV(:, :, 3) = RGB(:, :, 1).* A(3, 1) + RGB(:, :, 2).* A(3, 2) + RGB(:, :, 3).* A(3, 3);

%Add offset to U and V components, to avoid negative values
YUV(:, :, 2:3) = YUV(:, :, 2:3) + single(128); 

end

