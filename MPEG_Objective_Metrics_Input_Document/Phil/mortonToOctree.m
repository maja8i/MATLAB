function octreeBytes = mortonToOctree(mortonCodes,depth)
% MORTONTOOCTREE generates octree occupancy bytes from list of Morton codes
%
% mortonCodes = Nx1 (uint64) list of unique Morton codes, sorted
%
% depth = depth of octree, i.e., (number of bits in Morton code)/3
%
% octreeBytes = occupancy bytes in breadth-first order

% Generate tables by binary layer.
[I, W, F] = RAHTPrologue(mortonCodes,depth);

% Extract every third level of indexing, and add root.
III = vertcat({[1]},I(3*depth-2:-3:1)); % = cell(depth+1,1)

% Set flag to indicate if a code is a first child of its parent.
firstChildOfParentFlag = zeros(size(mortonCodes));

% Process one level at a time, from root (level=1) to leaves (level=depth+1)
octreeBytes = uint8.empty;
for level = 1:depth % level of parents
    
    % List Morton codes at child level.
    Mchildren = mortonCodes(III{level+1});
    
    % Extract 3-bit code (path to child from its parent).
    threeBitCode = bitand(bitshift(Mchildren,3*(level-depth)),uint64(7));
    
    % Find the (index of the) parent for each child.
    firstChildOfParentFlag(III{level}) = 1; % flag first child of each parent
    parentPtr = cumsum(firstChildOfParentFlag(III{level+1}));
    
    % Compute nodeByte for each parent.
    parentCount = length(III{level}); % count of parents
    childOccupiedFlag = zeros(8*parentCount,1,'uint8'); % flags if b'th child occupied
    childOccupiedFlag(uint64(8*(parentPtr-1))+threeBitCode+1) = 1; % 1 = child occupied
    childOccupiedFlag = reshape(childOccupiedFlag,8,parentCount); % indexed by (b,parent)
    nodeByte = zeros(parentCount,1,'uint8'); % allocate occupancy byte for each parent
    for b = 1:8 % child index
        nodeByte = bitor(nodeByte,bitshift(childOccupiedFlag(b,:)',8-b));
    end
    octreeBytes = [octreeBytes; nodeByte];
end
