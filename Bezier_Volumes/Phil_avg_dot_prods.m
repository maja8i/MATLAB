clear all
close all

depth = 10;
level = 7;

unitCube = [
    0 0 0;
    1 0 0;
    0 1 0;
    1 1 0;
    0 0 1;
    1 0 1;
    0 1 1;
    1 1 1];

% Create list of corners.
blockOrigin = [8 8 8; 16 16 16]; % (single(ot.SpatialIndex{level+1}*(2^(depth-level))) - 0.5);
blockOriginX8 = reshape(shiftdim(repmat(blockOrigin,1,1,8),2),8*size(blockOrigin,1),3);
blockCornerX8 = blockOriginX8 + repmat(unitCube*2^(depth-level),size(blockOrigin,1),1);

% Find delta between block corners and centroids.
blockCentroidX8 = blockOriginX8 + repmat([.5 .5 .5]*2^(depth-level),8*size(blockOrigin,1),1); % should use centroid, not corner
blockDeltaX8 = blockCornerX8 - blockCentroidX8;

% Get normals for each block and take dot product with deltas.
blockNormal = [1 -1 1; .9 -1 1.1]; % get block normals
blockNormalX8 = reshape(shiftdim(repmat(blockNormal,1,1,8),2),8*size(blockNormal,1),3);
dotProdX8 = sum(blockDeltaX8 .* blockNormalX8,2);

% Find list of unique corners, uniqueCorner, s.t. blockCornerX8 = uniqueCorner(iseg,:).
[uniqueCorner, ~, iCorner] = unique(blockCornerX8,'rows','stable');

% Create MxN matrix A, where M = length(uniqueCorner) and N = length(blockCornerX8),
% such that A(n,m) = 1 if blockCornerX8(n,:)==uniqueCorner(m,:) and = 0 otherwise.
M = size(uniqueCorner,1); % number of unique corners
N = size(blockCornerX8,1); % number of corners including repetitions from all blocks
A = sparse(iCorner,[1:N]',ones(N,1),M,N);

% Compute average dot product for each unique corner.
uniqueCount = full(sum(A,2));
uniqueTotal = A * dotProdX8;
uniqueAverage = single(uniqueTotal ./ repmat(uniqueCount,1,size(uniqueTotal,2))); % valid only when uniqueCount>0, else NaN