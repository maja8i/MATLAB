classdef octreeClass
    % OCTREECLASS represents an octree
    %   Detailed explanation goes here
    
    properties
        
        % Public
        Depth % depth of octree
        NodeCount % number of blocks at each level
        SpatialIndex % [x,y,z] index of block
        PathPrefix % interleaving of x, y, z block indices
        ParentPtr % pointer to parent
        FirstChildPtr % pointer to first child
        ChildCount % number of children
        FirstDescendantPtr % pointer to first descendant (occupied voxel) in block
        DescendantCount % number of descendants (occupied voxels) in block
        OccupancyCode % indicator of occupied children
        
        % Private
        HashTable
    end
    
    methods
        
        % Constructor.
        function obj = octreeClass(mortonCodes,depth)
            
            % Set defaults.
            obj.Depth = depth;
            obj.NodeCount = zeros(depth+1,1);
            obj.SpatialIndex = cell(depth+1,1);
            obj.PathPrefix = cell(depth+1,1);
            obj.ParentPtr = cell(depth+1,1);
            obj.FirstChildPtr = cell(depth,1);
            obj.ChildCount = cell(depth,1);
            obj.FirstDescendantPtr = cell(depth+1,1);
            obj.DescendantCount = cell(depth+1,1);
            obj.OccupancyCode = cell(depth,1);
            obj.HashTable = cell(depth+1,1);
            
            % Generate tables by binary layer.
            [I, W, F] = RAHTPrologue(mortonCodes,depth);

            % Extract every third level, and add root.
            obj.FirstDescendantPtr = vertcat({[1]},I(3*depth-2:-3:1));
            obj.DescendantCount = vertcat({[length(I{1})]},W(3*depth-2:-3:1));
            
            % Initialize root node.
            obj.NodeCount(1) = 1;
            obj.SpatialIndex{1} = [uint16(0),uint16(0),uint16(0)];
            obj.PathPrefix{1} = [uint64(0)];
            obj.ParentPtr{1} = [uint64(0)];
            obj.HashTable{1} = [uint64(1),uint64(0)];
            
            % Set flag to indicate if a code is a first child of its parent.
            firstChildOfParentFlag = zeros(size(mortonCodes));
            
            % Fill in tables by level.
            for level = 1:depth
                
                % Get Morton code for each child.
                M = mortonCodes(obj.FirstDescendantPtr{level+1});
                
                % Compute spatial index for each child.
                siX = zeros(length(M),1,'uint16');
                siY = zeros(length(M),1,'uint16');
                siZ = zeros(length(M),1,'uint16');
                for b = 1:level
                    siX = bitset(siX,b,bitget(M,3*(depth-level+b-1)+3));
                    siY = bitset(siY,b,bitget(M,3*(depth-level+b-1)+2));
                    siZ = bitset(siZ,b,bitget(M,3*(depth-level+b-1)+1));
                end
                obj.SpatialIndex{level+1} = [siX,siY,siZ];

                % Compute path prefix for each child.
                obj.PathPrefix{level+1} = bitshift(M,3*(level-depth));

                % Compute parent pointer for each child.
                firstChildOfParentFlag(obj.FirstDescendantPtr{level}) = 1; % flag first child of each parent
                obj.ParentPtr{level+1} = cumsum(firstChildOfParentFlag(obj.FirstDescendantPtr{level+1}));
                
                % Compute total number of children at each level.
                obj.NodeCount(level+1) = length(obj.ParentPtr{level+1});
                
                % Compute occupancy byte and child count for parent.
                % MSB indicates whether child with 3bitCode=000 is occupied or not.
                % LSB indicates whether child with 3bitCode=111 is occupied or not.
                threeBitCode = bitand(obj.PathPrefix{level+1},uint64(7)); % paths from parent to each child
                parentCount = length(obj.DescendantCount{level}); % count of parents (nodes at this level)
                childOccupiedFlag = zeros(8*parentCount,1,'uint8'); % flags if b'th child occupied
                childOccupiedFlag(uint64(8*(obj.ParentPtr{level+1}-1))+threeBitCode+1) = 1; % 1 = child occupied
                childOccupiedFlag = reshape(childOccupiedFlag,8,parentCount); % indexed by (b,parent)
                oB = zeros(parentCount,1,'uint8'); % allocate occupancy byte for each parent
                cC = zeros(parentCount,1,'uint8'); % allocate child count for each parent
                for b = 1:8 % child index
                    oB = bitor(oB,bitshift(childOccupiedFlag(b,:)',8-b));
                    cC = cC + childOccupiedFlag(b,:)';
                end
                obj.OccupancyCode{level} = oB;
                obj.ChildCount{level} = cC;
                
                % Compute pointer to first child.
                lastChildPtr = cumsum(int32(obj.ChildCount{level}));
                obj.FirstChildPtr{level} = uint32(lastChildPtr - int32(obj.ChildCount{level}) + 1);
                
                % Compute hash table.
                % First column is the key: prefix of the node in the octree
                % Second column is the value: pointer to the node in the level
                %hT = zeros(2^ceil(log2(3*length(M))),2,'uint64');
                hTsize = 2^ceil(log2(3*length(M))); % next power of 2 up from 3*length(M)
                hT = zeros(hTsize,2,'uint64');
                mask = uint64(hTsize - 1);
                numBits = level; % number of relevant bits in each spatial index
                key = xyzToMorton(obj.SpatialIndex{level+1},numBits);
                value = (1:length(M))';
                hash = bitand(key,mask) + 1; % simple hash
                %hash = bitand(bitxor(key,key.*key),mask) + 1; % better hash
                %hash = uint64(bitxor(key,bitshift(key,-33)) * uint64(18397678908172766413));
                %hash = uint64(bitxor(hash,bitshift(hash,-33)) * uint64(14181476777654086739));
                %hash = uint64(bitand(bitxor(hash,bitshift(hash,-33)),mask) + 1);
%                 meanBucketSize = length(M) / sum(unique(hash)~=0);
%                 meanProbeCount = 0;
%                 sparseTable = sparse(value,double(hash),ones(length(M),1),length(M),hTsize);
%                 maxBucketSize = max(ones(1,length(M))*sparseTable);
                
                probeCount = 0;
                while ~isempty(key)
                    probeCount = probeCount + 1;
                    
                    % Sort on hash, so it's easier to determine duplicates
                    [hash, I] = sort(hash);
                    key = key(I);
                    value = value(I);
                    
                    % No collision as long as it's hash is not equal to the previous hash
                    % and its entry in the hash table is not already filled.
                    noCollision = (hash ~= [0;hash(1:end-1)]) & (hT(hash,2) == 0);
%                     meanProbeCount = meanProbeCount + probeCount * sum(noCollision);
                    
                    % Put those without collisions in the hash table.
                    hT(hash(noCollision),:) = [key(noCollision), value(noCollision)];
                    
                    % Winnow the list down to the collisions, and try again.
                    collision = ~noCollision;
                    key = key(collision);
                    value = value(collision);
                    %hash = bitand(hash(collision),mask) + 1; % linear probing
                    hash = bitand(hash(collision)+probeCount-1,mask) + 1; % quadratic probing
                end
%                 meanProbeCount = meanProbeCount / length(M);
%                 fprintf('level=%d meanBucketSize=%g maxBucketSize=%d meanProbeCount=%g maxProbeCount=%d\n',level,meanBucketSize,maxBucketSize,meanProbeCount,probeCount);
                
                obj.HashTable{level+1} = hT;
                
            end
        end
        
        % Return 'depth' levels of occupancy codes concatenated in breadth-first order.
        function oC = occupancyCodes(obj,depth)
            if nargin < 2
                depth = obj.Depth;
            end
            oC = uint8.empty;
            for level = 1:depth
                oC = [oC; obj.OccupancyCode{level}];
            end
        end
        
        % Return pointer to node having SpatialIndex.
        % If no such node exists in the level, return 0.
        % SpatialIndex is a Nx3 uint16 array.
        % 'depth' is number of relevant bits in each SpatialIndex component.
        % Pointers are returned in a Nx1 uint64 array.
        function nP = nodePtr(obj,spatialIndex,depth)
            mask = uint64(size(obj.HashTable{depth+1},1)-1);
            key = xyzToMorton(spatialIndex,depth);
            hash = bitand(key,mask) + 1; % simple hash
            %hash = bitand(bitxor(key,key.*key),mask) + 1; % better hash
            key_value = obj.HashTable{depth+1}(hash,1:2);
            collision = (key_value(:,2) ~= 0) & (key_value(:,1) ~= key);
            probeCount = 0;
            while any(collision)
                probeCount = probeCount + 1;
                %hash(collision) = bitand(hash(collision),mask) + 1; % linear probing
                hash(collision) = bitand(hash(collision)+probeCount-1,mask) + 1; % quadratic probing
                key_value(collision,:) = obj.HashTable{depth+1}(hash(collision),1:2);
                collision = (key_value(:,2) ~= 0) & (key_value(:,1) ~= key);
            end
            nP = key_value(:,2);
        end
    end
    
end

