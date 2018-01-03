function [transformedAttributes, weights] = RAHT(I, W, F, attributes, depth, fracStepsize, bHaarThresh)
%RAHT: Region Adaptive Hierarchical Transform
%
% Copyright 8i Labs, Inc., 2017
% This code is to be used solely for the purpose of developing the MPEG PCC standard.
%
%   Transform attributes into transformedAttributes.
%
%   Input:
%   I, W, F are all 3*depth x 1 cell arrays
%   I{b} = indices of Morton codes of lead voxels at binary level b
%   W{b} = weights of cells at binary level b
%   F{b} = boolean flags to indicate left siblings at binary level b
%   attributes = NxK array of attributes
%   depth = number of bits for each of x, y, and z
%   fracStepsize = stepsize for quantizing fraction of left descendants
%       or zero if no quantization
%   bHaarThresh = threshold for b, below which transform is Haar
%
%   Output:
%   transformedAttributes = NxK array of transformed attributes
%   weights = Nx1 array of weights of transformed attributes
%
%   Note:
%   level b=1 is low order bit of Morton code, or leaves of binary tree
%   level b=3*depth is high order bit of Morton code, or root of tree

% Set defaults.
if nargin < 6
    fracStepsize = 0; % no quantization
end
if nargin < 7
    bHaarThresh = 0; % no threshold
end

% Initialize transformed attributed to attributes,
% because we're going to transform in place.
transformedAttributes = attributes;

% Initialize weights to all ones.
weights = ones(size(transformedAttributes,1),1);

% Process one level at a time, from leaves (b=1) to root (b=3*depth).
for b = 1:3*depth
    i0 = I{b}([F{b};0] == 1); % left sibling indices
    if isempty(i0)
        continue;
    end
    i1 = I{b}([0;F{b}] == 1); % right sibling indices
    w0 = W{b}([F{b};0] == 1); % left sibling weights
    w1 = W{b}([0;F{b}] == 1); % right sibling weights
    x0 = transformedAttributes(i0,:); % left sibling coefficients
    x1 = transformedAttributes(i1,:); % right sibling coefficients
    if fracStepsize == 0
        frac = w0 ./ (w0+w1);
    else
        frac = round((w0 ./ (w0+w1))/fracStepsize) * fracStepsize;
    end
    if b <= bHaarThresh
        frac(frac>0 & frac<1) = 0.5;
    end
    alpha = repmat(single(sqrt(frac)),1,size(transformedAttributes,2));
    beta = repmat(single(sqrt(1-frac)),1,size(transformedAttributes,2));
    transformedAttributes(i0,:) = alpha .* x0 + beta .* x1;
    transformedAttributes(i1,:) = -beta .* x0 + alpha .* x1;
    weights(i0) = weights(i0) + weights(i1);
    weights(i1) = weights(i0); % same as i0
%     var0 = weights(i0);
%     var1 = weights(i1);
%     alpha1 = alpha(:,1);
%     beta1 = beta(:,1);
%     crossterm = 2*((0.95)^(2^floor((b-1)/3)))*alpha1.*beta1.*sqrt(var0.*var1);
%     weights(i0) = (alpha1.^2).*var0 + crossterm + (beta1.^2).*var1;
%     weights(i1) = (beta1.^2).*var0 - crossterm + (alpha1.^2).*var1;
end
