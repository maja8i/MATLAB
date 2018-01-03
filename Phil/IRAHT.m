function attributes = IRAHT(I, W, F, transformedAttributes, depth, level, fracStepsize, bHaarThresh)
%RAHT: Inverse Region Adaptive Hierarchical Transform
%
% Copyright 8i Labs, Inc., 2017
% This code is to be used solely for the purpose of developing the MPEG PCC standard.
%
%   Inverse transform transformedAttributes into attributes.
%
%   Input:
%   I, W, F are all 3*depth x 1 cell arrays
%   I{b} = indices of Morton codes of lead voxels at binary level b
%   W{b} = weights of cells at binary level b
%   F{b} = boolean flags to indicate left siblings at binary level b
%   transformedAttributes = NxK array of transfomed attributes
%   depth = number of bits for each of x, y, and z
%   level (optional) = level to which the inverse transform is performed
%       default is level = depth
%   fracStepsize = stepsize for quantizing fraction of left descendants
%       or zero if no quantization
%   bHaarThresh = threshold for b, below which transform is Haar
%
%   Output:
%   attributes = NxK array of attributes
%
%   Note:
%   level b=1 is low order bit of Morton code, or leaves of binary tree
%   level b=3*depth is high order bit of Morton code, or root of tree

% Set defaults.
if nargin < 6
    level = depth;
end
if nargin < 7
    fracStepsize = 0; % no quantization
end
if nargin < 8
    bHaarThresh = 0; % no threshold
end

% Initialize attributed to transformed attributes,
% because we're going to inverse transform in place.
attributes = transformedAttributes;

% Process one level at a time, from root (b=3*depth) to leaves (b=1).
for b = 3*depth:-1:(3*(depth-level)+1)
    i0 = I{b}([F{b};0] == 1); % left sibling indices
    if isempty(i0)
        continue;
    end
    i1 = I{b}([0;F{b}] == 1); % right sibling indices
    w0 = W{b}([F{b};0] == 1); % left sibling weights
    w1 = W{b}([0;F{b}] == 1); % right sibling weights
    x0 = attributes(i0,:); % left sibling coefficients
    x1 = attributes(i1,:); % right sibling coefficients
    if fracStepsize == 0
        frac = w0 ./ (w0+w1);
    else
        frac = round((w0 ./ (w0+w1))/fracStepsize) * fracStepsize;
    end
    if b <= bHaarThresh
        frac(frac>0 & frac<1) = 0.5;
    end
    alpha = repmat(single(sqrt(frac)),1,size(attributes,2));
    beta = repmat(single(sqrt(1-frac)),1,size(attributes,2));
    attributes(i0,:) = alpha .* x0 - beta .* x1;
    attributes(i1,:) = beta .* x0 + alpha .* x1;
end