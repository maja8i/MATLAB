% Read voxelizedPointCloud from file.
%
function obj = read(obj,filename)

    fprintf('reading filename=%s...\n',filename);
    
    % Process by file type.
    [pathstr,name,ext] = fileparts(filename);
    switch lower(ext)
        case '.ply'
            % Read into plyStruct.
            ply = plyRead(filename);

            % Pay attention only to the vertex element.
            i = find(strcmpi('vertex',ply.elemNameList));

            % Get number of properties per vertex.
            propCount = length(ply.propNameListList{i});

            % Traverse properties and fill in slots.
            % Slots are partitioned into location, color, and normal slots.
            % Each of these is further partitioned into 1st, 2nd, and 3rd components.
            % (Color may have a 4th component.)
            % For location, components are x, y, z.
            % For color, components are red, green, blue; Y, U, V; or Y, Cb, Cr.
            % (A fourth color component may be alpha.)
            % For normal, components are nx, ny, nz.
            % Color slots may also have a camera direction given by a number in the name.
            % For example, red1, ..., red14 are the red components for cameras 1, ..., 14.
            attributeCount = 0; % number of generic attributes
            for j = 1:propCount

                % Fill in the slot according to name.
                name = ply.propNameListList{i}(j);
                switch name
                    case 'x' %
                        obj.Location = eval([class(ply.propArrayListList{i}{j}) '(obj.Location)']);
                        obj.Location(:,1) = ply.propArrayListList{i}{j};
                    case 'y'
                        obj.Location = eval([class(ply.propArrayListList{i}{j}) '(obj.Location)']);
                        obj.Location(:,2) = ply.propArrayListList{i}{j};
                    case 'z'
                        obj.Location = eval([class(ply.propArrayListList{i}{j}) '(obj.Location)']);
                        obj.Location(:,3) = ply.propArrayListList{i}{j};
                    case 'nx'
                        obj.Normal = eval([class(ply.propArrayListList{i}{j}) '(obj.Normal)']);
                        obj.Normal(:,1) = ply.propArrayListList{i}{j};
                    case 'ny'
                        obj.Normal = eval([class(ply.propArrayListList{i}{j}) '(obj.Normal)']);
                        obj.Normal(:,2) = ply.propArrayListList{i}{j};
                    case 'nz'
                        obj.Normal = eval([class(ply.propArrayListList{i}{j}) '(obj.Normal)']);
                        obj.Normal(:,3) = ply.propArrayListList{i}{j};
                    case {'red', 'R', 'Y'}
                        obj.Color(:,1) = single(ply.propArrayListList{i}{j}); % force to single
                    case {'green', 'G', 'U', 'Cb'}
                        obj.Color(:,2) = single(ply.propArrayListList{i}{j}); % force to single
                    case {'blue', 'B', 'V', 'Cr'}
                        obj.Color(:,3) = single(ply.propArrayListList{i}{j}); % force to single
                    case {'alpha', 'A'}
                        obj.Color(:,4) = single(ply.propArrayListList{i}{j}); % force to single
                    otherwise
                        % Maybe this is a view-dependendent color
                        % and there is a camera number in the name.
                        % Try to split it out.
                        [nameCell, numberCell] = regexp(name,'\d*','split','match');
                        if ~isempty(nameCell) && ~isempty(numberCell)
                            n = str2double(numberCell{1});
                            switch nameCell{1}
                                case {'red', 'Y'}
                                    obj.ViewDependentColor(:,3*(n-1)+1) = single(ply.propArrayListList{i}{j});
                                case {'green', 'U', 'Cb'}
                                    obj.ViewDependentColor(:,3*(n-1)+2) = single(ply.propArrayListList{i}{j});
                                case {'blue', 'V', 'Cr'}
                                    obj.ViewDependentColor(:,3*(n-1)+3) = single(ply.propArrayListList{i}{j});
                            end
                        else % This is a generic attribute.
                            attributeCount = attributeCount + 1;
                            obj.Attributes(:,attributeCount) = single(ply.propArrayListList{i}{j}); % force to single
                            obj.AttributeNames(attributeCount) = name;
                            obj.AttributeTypes(attributeCount) = class(ply.propArrayListList{i}{j}); % keep track of type
                        end
                end

                % Set ColorSpace.
                switch name{1:end}
                    case {'red', 'green', 'blue'}
                        obj.ColorSpace = 'RGB';
                    case {'U', 'V'}
                        obj.ColorSpace = 'YUV';
                    case {'Cb', 'Cr'}
                        obj.ColorSpace = 'YCbCr';
                end
            end

            % Set Count and update Limits if necessary.
            obj = updateCountAndLimits(obj);
            
            % Parse ply comments, if any.
            for i = 1:length(ply.comments)
                comment = lower(ply.comments(i));
                s = split(comment);
                if length(s) >= 2
                    switch s(2)
                        case 'frame_to_world_scale'
                            if length(s) >= 3
                                obj.FrameToWorldScale = str2double(s(3));
                            end
                        case 'frame_to_world_translation'
                            if length(s) >= 5
                                obj.FrameToWorldTranslation = str2double(s(3:5));
                            end
                        case 'width'
                            if length(s) >= 3
                                obj.CubeWidth = str2double(s(3));
                            end
                    end
                end
            end
        otherwise
            error('cannot read filetype=%s\n',ext);
    end
end