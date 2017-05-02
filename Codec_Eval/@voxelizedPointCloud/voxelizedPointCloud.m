% pointCloud class
classdef voxelizedPointCloud % < pointCloud
    properties
        % properties of standard matlab pointCloud
        Location
        Normal
        Color
        Count
        XLimits
        YLimits
        ZLimits
        % extensions
        ViewDependentColor % Count-by-(3*ViewCount) array
        Attributes % Count-by-attributeCount array of attributes
        AttributeNames % 1-by-attributeCount string array of attribute names
        AttributeTypes % 1-by-attributeCount string array of attribute names
        ColorSpace % 'RGB', 'YUV', 'YCbCr'
        FrameToWorldScale % scalar transformation from frame to world coords
        FrameToWorldTranslation % vector translation from frame to world coords
        CubeWidth % frame coords lie in [0,CubeWidth]
        Depth % cube is voxelized into 2^depth voxels per side
        Morton % list of Morton codes
%         TransformedColor % transformed color
%         TransformedColorWeights % transformed color weights
%         TransformedColorIndices % indices of quantized transformed colors
%         TransformedColorReproductions % quantized transformed colors
    end
    methods
        
        % Constructor.
        function obj = voxelizedPointCloud(filename)
            
            % Set defaults.
            obj.Location = single.empty;
            obj.Normal = single.empty;
            obj.Color = single.empty;
            obj.Count = 0;
            obj.ViewDependentColor = single.empty;
            obj.Attributes = single.empty;
            obj.AttributeNames = string.empty;
            obj.AttributeTypes = string.empty;
            obj.XLimits = [];
            obj.YLimits = [];
            obj.ZLimits = [];
            obj.ColorSpace = '';
            obj.FrameToWorldScale = 1.0;
            obj.FrameToWorldTranslation = [0 0 0];
            obj.CubeWidth = Inf;
            obj.Depth = 0;
            obj.Morton = uint64.empty;
%             obj.TransformedColor = single.empty;
%             obj.TransformedColorWeights = single.empty;
%             obj.TransformedColorIndices = single.empty;
%             obj.TransformedColorReproductions = single.empty;
            
            if nargin > 0
                obj = read(obj,filename);
            end
            
        end
        
        % Read from file.
        obj = read(obj,filename)
        
        % Write to file.
        write(obj,filename,format)
        
        % Remove all entries for which x, y, or z equal NaN or Inf.
        function obj = removeInvalidPoints(obj)
            if ~isempty(obj.Location)
                keep = (sum(isfinite(obj.Location),2) == 3);
                obj.Location = obj.Location(keep,:);
                if ~isempty(obj.Normal)
                    obj.Normal = obj.Normal(keep,:);
                end
                if ~isempty(obj.Color)
                    obj.Color = obj.Color(keep,:);
                end
                if ~isempty(obj.ViewDependentColor)
                    obj.ViewDependentColor = obj.ViewDependentColor(keep,:);
                end
                if ~isempty(obj.Attributes)
                    obj.Attributes = obj.Attributes(keep,:);
                end
                obj = updateCountAndLimits(obj);
            end
        end
        
        % Update count and limits.
        function obj = updateCountAndLimits(obj)
            if ~isempty(obj.Location)
                obj.Count = length(obj.Location(:,1));
                minLocation = min(obj.Location);
                maxLocation = max(obj.Location);
                obj.XLimits = [minLocation(1) maxLocation(1)];
                obj.YLimits = [minLocation(2) maxLocation(2)];
                obj.ZLimits = [minLocation(3) maxLocation(3)];
            end
        end
        
        % Set FrameToWorld transform
        % given 'cubeWidth', 'lowerLim', and 'maxWidth'
        % so that the FrameToWorld transform maps the frame-coordinate cube
        %    [0,cubeWidth]^3
        % to the world-coordinate cube
        %    [lowerLim(1),lowerLim(1)+maxWidth]
        %  x [lowerLim(2),lowerLim(2)+maxWidth]
        %  x [lowerLim(3),lowerLim(3)+maxWidth].
        % Argument 'cubeWidth' is required but arbitrary, depending on the
        % desired range of the frame coords, which should lie in [0,cubeWidth]^3.
        % ('cubeWidth' is usually either 1.0 or 2^depth-1, depth is octree depth.)
        % Arguments 'lowerLim' and 'maxWidth' are optional, but if specified,
        % should be specified such that all world coords are contained the cube
        %    [lowerLim(1),lowerLim(1)+maxWidth]
        %  x [lowerLim(2),lowerLim(2)+maxWidth]
        %  x [lowerLim(3),lowerLim(3)+maxWidth].
        % If not specified, 'lowerLim' and 'maxWidth' are derived from
        % the minimal enclosing bounding cube of the world coordinates.
        function obj = setTransform(obj,cubeWidth,lowerLim,maxWidth)
            if nargin <= 2 % lowerLim and maxWidth not specified
                % Derive from minimal enclosing bounding cube.
                lowerLim = [obj.XLimits(1) obj.YLimits(1) obj.ZLimits(1)];
                upperLim = [obj.XLimits(2) obj.YLimits(2) obj.ZLimits(2)];
                maxWidth = max(upperLim - lowerLim);
            end
            % Set FrameToWorld transform.
            obj.FrameToWorldScale = maxWidth / cubeWidth;
            obj.FrameToWorldTranslation = lowerLim;
            obj.CubeWidth = cubeWidth;
        end
        
        % Frame coords nominally lie in [0,CubeWidth]^3,
        % while world coords can be anything.
        
        % Transform (x,y,z) from world coords to frame coords.
        function obj = worldToFrame(obj)
            obj.Location = (obj.Location - repmat(obj.FrameToWorldTranslation,obj.Count,1)) ...
                / obj.FrameToWorldScale;
            obj = updateCountAndLimits(obj);
        end
        
        % Transform (x,y,z) from frame coords to world coords.
        function obj = frameToWorld(obj)
            obj.Location = obj.Location * obj.FrameToWorldScale ...
                + repmat(obj.FrameToWorldTranslation,obj.Count,1);
            obj = updateCountAndLimits(obj);
        end
        
        % Compute Morton codes by integerizing x, y, and z and interleaving.
        % Then sort.
        function obj = mortonizeAndSort(obj)
            
            % Allocate space for Morton codes.
            obj.Morton = zeros(obj.Count,1,'uint64');
            
            % Integerize x, y, z Locations.
            intLocation = uint64(obj.Location);
            
            % Multiplex bits of integerized x, y, z into Morton codes.
            for b = 1:21
                obj.Morton = bitor(obj.Morton,bitshift(bitget(intLocation(:,1),b),(b-1)*3+2));
                obj.Morton = bitor(obj.Morton,bitshift(bitget(intLocation(:,2),b),(b-1)*3+1));
                obj.Morton = bitor(obj.Morton,bitshift(bitget(intLocation(:,3),b),(b-1)*3));
            end
            
            % Sort Morton codes.
            [obj.Morton, I] = sort(obj.Morton);
            for k = 1:3
                obj.Location(:,k) = obj.Location(I,k);
            end
            if ~isempty(obj.Normal)
                for k = 1:3
                    obj.Normal(:,k) = obj.Normal(I,k);
                end
            end
            if ~isempty(obj.Color)
                for k = 1:size(obj.Color,2)
                    obj.Color(:,k) = obj.Color(I,k);
                end
            end
            if ~isempty(obj.ViewDependentColor)
                for k = 1:size(obj.ViewDependentColor,2)
                    obj.ViewDependentColor(:,k) = obj.ViewDependentColor(I,k);
                end
            end
            if ~isempty(obj.Attributes)
                for k = 1:size(obj.Attributes,2)
                    obj.Attributes(:,k) = obj.Attributes(I,k);
                end
            end
        end
        
        % Voxelize by averaging all values within the same Morton Code.
        function obj = voxelize(obj)
            
            % Make sure Morton codes have been created and sorted.
            if isempty(obj.Morton)
                obj = mortonizeAndSort(obj);
            end
            
            % Indicate 1 if Morton code is different from previous Morton code in list.
            differentFromPrevious = [1; obj.Morton(2:end) ~= obj.Morton(1:end-1)];
            
            % Find list of unique Morton codes, uniqueMorton,
            % i.e., [uniqueMorton, iMorton, ~] = unique(Morton),
            % s.t. uniqueMorton = Mortion(iMorton).  Equivalently:
            uniqueMorton = obj.Morton(logical(differentFromPrevious));
            iMorton = cumsum(differentFromPrevious);
            
            % Create MxN matrix A, where M = length(Morton) and N = length(uniqueMorton),
            % such that A(m,n) = 1 if Morton(m)==uniqueMorton(n) and = 0 otherwise.
            M = length(obj.Morton); % same as obj.Count
            N = iMorton(end); % same as length(uniqueMorton)
            A = sparse([1:M]',iMorton,ones(M,1),M,N);
            
            % Average Normal, Color, and ViewDependentColor for each uniqueMorton code.
            uniqueCounts = (ones(1,M) * A)'; % how many times each uniqueMorton code is used
            if ~isempty(obj.Normal)
                obj.Normal = (A' * obj.Normal) ./ uniqueCounts;
            end
            if ~isempty(obj.Color)
                obj.Color = single((A' * double(obj.Color)) ./ uniqueCounts);
            end
            if ~isempty(obj.ViewDependentColor)
                obj.ViewDependentColor = (A' * obj.ViewDependentColor) ./ uniqueCounts;
            end
            if ~isempty(obj.Attributes)
                obj.Attributes = (A' * double(obj.Attributes)) ./ uniqueCounts;
            end
            obj.Morton = uniqueMorton;
            obj = deMortonize(obj);
        end
        
        % De-mortonize.
        function obj = deMortonize(obj)
            obj.Count = length(obj.Morton);
            intLocation = zeros(obj.Count,3,'uint64');
            
            % Demultiplex bits of Morton codes in to x, y, z Locations.
            for b = 0:20
                intLocation(:,1) = bitor(intLocation(:,1),bitshift(bitget(obj.Morton,b*3+3),b));
                intLocation(:,2) = bitor(intLocation(:,2),bitshift(bitget(obj.Morton,b*3+2),b));
                intLocation(:,3) = bitor(intLocation(:,3),bitshift(bitget(obj.Morton,b*3+1),b));
            end
            obj.Location = single(intLocation);
        end
        
        % Transform.
        function [transformedAttributes, weights] = transform(obj,attributes)
            
            % Make sure Morton codes have been created and sorted.
            if isempty(obj.Morton)
                error('mortonizeAndSort not yet called');
            end
            
            [I, W, F] = RAHTPrologue(obj.Morton,obj.Depth);
            [transformedAttributes, weights] = RAHT(I,W,F,attributes,obj.Depth);

        end
        
        % Inverse transform.
        function attributes = invTransform(obj,transformedAttributes,level)
            
            % Level to which we inverse transform is optional.
            if nargin < 3
                level = obj.Depth;
            end
            
            % Make sure Morton Codes have been created and sorted.
            if isempty(obj.Morton)
                error('mortonizeAndSort not yet called');
            end
            
            [I, W, F] = RAHTPrologue(obj.Morton,obj.Depth);
            attributes = IRAHT(I,W,F,transformedAttributes,obj.Depth,level);

        end
        
        % Project points orthogonally onto image.
        % Degrees is amount to rotate object clockwise looking down along Y axis.
        function cdata = orthogonalProjection(obj,degrees)
            
            % Points lie in [0,cubeWidth]^3.
            cubeWidth = obj.CubeWidth;
            
            % Get points, but flip coordinates as necessary.
            X = obj.Location(:,1);
            Y = cubeWidth-obj.Location(:,2); % flip Y for image
            Z = cubeWidth-obj.Location(:,3); % flip Z for depthmap

            % Rotate around Y axis if necessary.
            XZ = [X, Z] - cubeWidth/2;
            theta = degrees*(pi/180);
            R = [cos(theta), -sin(theta); sin(theta),  cos(theta)]; % rotation matrix
            XZ = (XZ * R) + cubeWidth/2;
            XZ = round(XZ);
            X = XZ(:,1);
            Z = XZ(:,2);

            % Implement depth buffer, keeping largest Z for same X, Y.
            A = sortrows([X, Y, Z, obj.Color]);
            I = ((A(:,1)~=[0;A(1:end-1,1)]) | (A(:,2)~=[0;A(1:end-1,2)])) & A(:,1)>=0 & A(:,1)<=cubeWidth;
            A = A(I,:);

            % Fill image.
            XY = A(:,1) * (cubeWidth+1) + A(:,2) + 1;
            cdata = zeros((cubeWidth+1)*(cubeWidth+1),3,'uint8');
            cdata(XY,:) = uint8(A(:,4:6));
            cdata = reshape(cdata,(cubeWidth+1),(cubeWidth+1),3);
        end
    end
end