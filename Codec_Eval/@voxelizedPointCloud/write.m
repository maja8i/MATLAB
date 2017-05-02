% pcwrite
function write(obj,filename,format)

    if nargin < 3
        format = 'ascii';
    end
    
    fprintf('writing filename=%s format=%s...\n',filename,format);
    
    % Process by file type.
    [pathstr,name,ext] = fileparts(filename);
    switch lower(ext)
        case '.ply'
            
            % Construct plyStruct with single element.
            elemNameList = [string('vertex')];
            propNameList = string.empty;
            propTypeList = string.empty;
            propArrayList = cell.empty;
            comments = [string(sprintf('comment frame_to_world_scale %g',obj.FrameToWorldScale)); ...
                string(sprintf('comment frame_to_world_translation %g %g %g',obj.FrameToWorldTranslation));
                string(sprintf('comment width %g',obj.CubeWidth))];
            
            % Fill properties.
            if ~isempty(obj.Location)
                propNameList = [propNameList string('x') string('y') string('z')];
                propArrayList = [propArrayList obj.Location(:,1) obj.Location(:,2) obj.Location(:,3)];
                type = string(class(obj.Location));
                propTypeList = [propTypeList type type type];
            end
            if ~isempty(obj.Normal)                
                propNameList = [propNameList string('nx') string('ny') string('nz')];
                propArrayList = [propArrayList obj.Normal(:,1) obj.Normal(:,2) obj.Normal(:,3)];
                type = string(class(obj.Normal));
                propTypeList = [propTypeList type type type];
            end
            switch obj.ColorSpace
                case 'RGB'
                    col = {'red', 'green', 'blue', 'alpha'};
                case 'YUV'
                    col = {'Y', 'U', 'V', 'alpha'};
                case 'YCbCr'
                    col = {'Y', 'Cb', 'Cr', 'alpha'};
            end
            if ~isempty(obj.Color)
                for j = 1:size(obj.Color,2)
                    propNameList = [propNameList string(col{j})];
                    propArrayList = [propArrayList uint8(obj.Color(:,j))];
                    propTypeList = [propTypeList string('uint8')]; % force to uint8
                end
            end
            if ~isempty(obj.ViewDependentColor)
                for j = 1:size(obj.ViewDependentColor,2)
                    propNameList = [propNameList strcat(col{mod(j-1,3)+1},string(j))];
                    propArrayList = [propArrayList uint8(obj.ViewDependentColor(:,j))];
                    propTypeList = [propTypeList string('uint8')]; % force to uint8
                end
            end
            if ~isempty(obj.Attributes)
                for j = 1:size(obj.Attributes,2)
                    propNameList = [propNameList obj.AttributeNames(j)];
                    propArrayList = [propArrayList obj.Attributes(:,j)];
                    propTypeList = [propTypeList obj.AttributeTypes(j)];
                end
            end
            
            % Convert matlab types to ply types.
            for j = 1:length(propTypeList)
                switch propTypeList(j)
                    case 'int8'
                        propTypeList(j) = string('char');
                    case 'uint8'
                        propTypeList(j) = string('uchar');
                    case 'int16'
                        propTypeList(j) = string('short');
                    case 'uint16'
                        propTypeList(j) = string('ushort');
                    case 'int32'
                        propTypeList(j) = string('int');
                    case 'uint32'
                        propTypeList(j) = string('uint');
                    case 'single'
                        propTypeList(j) = string('float');
                    case 'double'
                        propTypeList(j) = string('double');
                end
            end
            
            propNameListList = {{propNameList}};
            propTypeListList = {{propTypeList}};
            propArrayListList = {{propArrayList}};
            plyStruct = struct('elemNameList',elemNameList,...
                'propNameListList',propNameListList,...
                'propTypeListList',propTypeListList,...
                'propArrayListList',propArrayListList,...
                'comments',comments);
            plyWrite(plyStruct,filename,format);
        otherwise
            error('cannot write filetype=%s\n',ext);            
    end
end