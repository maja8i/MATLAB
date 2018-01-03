function [I, W, F] = RAHTPrologue(M1,depth)
%RAHTPROLOGUE Prologue for RAHT and IRAHT
%
% Copyright 8i Labs, Inc., 2017
% This code is to be used solely for the purpose of developing the MPEG PCC standard.
%
% Prologue to RAHT and IRAHT
%
%   Input:
%   M1 = list of Morton codes, sorted
%   depth = number of bits for each of x, y, and z
%
%   Output:
%   I, W, F are all 3*depth x 1 cell arrays
%   I{b} = indices of Morton codes of lead voxels at binary level b
%   W{b} = weights of cells at binary level b
%   F{b} = boolean flags to indicate left siblings at binary level b
%
%   Note:
%   level b=1 is low order bit of Morton code, or leaves of binary tree
%   level b=3*depth is high order bit of Morton code, or root of tree

% Get number of voxels, or leaves of binary tree.
N = length(M1);

% Allocate cell arrays.
I = cell(3*depth,1);
W = cell(3*depth,1);
F = cell(3*depth,1);

c1 = ones(1,'uint64'); % constant 1

% Process one level at a time, from leaves (b=1) to root (b=3*depth).
for b = 1:3*depth
    if b == 1 % initialize indices of coeffs at level b
        I{b} = uint32(1:N)'; % vector of indices from 1 to N
    else % define indices of coeffs at level b
        I{b} = I{b-1}(~[0;F{b-1}]);
    end
    Mb = M1(I{b}); % Morton codes at level b
    W{b} = double([I{b}(2:end);N+1] - I{b}); % weights
    D = bitxor(Mb(1:end-1),Mb(2:end)); % path diffs
    F{b} = (bitand(D,(bitshift(c1,3*depth)-bitshift(c1,b)))) == 0; % is left sibling
end