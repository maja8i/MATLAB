function RGB = yuvToRgbImage(YUV)
% Convert YUV to RGB
%
% Converts YUV values at full scale to RGB values in [0,255]^3.
%   Thus Y takes values in the full range [0,255]
%   and U, V may take values outside this range.
%
% YUV = MxNx3 input array (single), where MxN are the dimensions of the 
%       input image
%
% RGB = MxNx3 output array (single)

YUV = single(YUV); % convert to floating point if necessary

A = single([    1,          0,          1.13983;...
                1,          -0.39465,   -0.58060;...
                1,          2.03211,    0]);

%Remove offset that was added to U and V components in the RGB -> YUV 
%conversion
YUV(:, :, 2:3) = YUV(:, :, 2:3) - single(128); 
            
RGB(:, :, 1) = YUV(:, :, 1).* A(1, 1) + YUV(:, :, 2).* A(1, 2) + YUV(:, :, 3).* A(1, 3);
RGB(:, :, 2) = YUV(:, :, 1).* A(2, 1) + YUV(:, :, 2).* A(2, 2) + YUV(:, :, 3).* A(2, 3);
RGB(:, :, 3) = YUV(:, :, 1).* A(3, 1) + YUV(:, :, 2).* A(3, 2) + YUV(:, :, 3).* A(3, 3);

end

