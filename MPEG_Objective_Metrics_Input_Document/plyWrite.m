% Write generic PLY file
% (arbitrary element and property names, arbitrary types; list types not implemented)
function plyWrite(ply,filename,format)
%
% ply.elemNameList - elemCount x 1 string array of element names,
%   e.g., ["vertex"; "face"], where "vertex" means string('vertex')
% ply.propNameListList - elemCount x 1 cell array of string arrays of property names,
%   e.g., {["x","y","z","red","green","blue"];["vertex_indices"]}
% ply.propTypeListList - elemCount x 1 cell array of string arrays of property types,
%   e.g., {["float","float","float","uchar","uchar","uchar"};["list uchar int"]}
% ply.propArrayListList - elemCount x 1 cell array of cell arrays of property arrays,
%   e.g., {{x,y,z,red,green,blue};{vertex_indices}}
%   where x, y, and z are nx1 float; red, green, and blue are nx1 uchar;
%   and vertex_indices are mxf int, where n and m are vertex and face counts
% ply.comments - k x 1 string array of comment lines
% Note 'list' types are not currently supported
%
% filename - character string including .ply extension
%
% format - 'ascii', 'binary_little_endian', or 'binary_big_endian'

% Typical header:
% ply
% format ascii 1.0 | format binary_big_endian 1.0
% element vertex 3073092
% property float x
% property float y
% property float z
% property uchar red
% property uchar green
% property uchar blue
% element face 0
% property list uchar int vertex_indices
% end_header

% Check format.
switch format
    case 'ascii'
    case 'binary_little_endian'
    case 'binary_big_endian'
    otherwise
        error('format "%s" not supported\n',format);
end

% Get number of elements.
elemCount = length(ply.elemNameList);

% Write file header.
fid = fopen(filename,'w');
if fid < 0
    error('cannot open file "%s" for writing\n',filename);
end
fprintf(fid,'ply\n');
fprintf(fid,'format %s 1.0\n',format);
for i = 1:length(ply.comments)
    fprintf(fid,'%s\n',ply.comments(i));
end
for i = 1:elemCount
    dataCount = size(ply.propArrayListList{i}{1},1);
    propCount = length(ply.propNameListList{i});
    fprintf(fid,'element %s %d\n',ply.elemNameList(i),dataCount);
    for j = 1:propCount
        fprintf(fid,'property %s %s\n',ply.propTypeListList{i}(j),ply.propNameListList{i}(j));
    end
end
fprintf(fid,'end_header\n');
headerByteCount = ftell(fid);
fclose(fid);

% Re-open file with the correct format.
switch format
    case 'ascii'
        fid = fopen(filename,'r+','n','US-ASCII');
    case 'binary_little_endian'
        fid = fopen(filename,'r+','ieee-le');
    case 'binary_big_endian'
        fid = fopen(filename,'r+','ieee-be');
    otherwise
        error('format "%s" not supported\n',format);
end

% Seek to beginning of data.
if fseek(fid,headerByteCount,'bof') ~= 0
    error('incomplete ply file');
end

% Write data for each element in turn.
elemFilePos = headerByteCount; % position of current element's data in file
for i = 1:elemCount
    propNameList = ply.propNameListList{i}; % list of property names for this element
    propTypeList = ply.propTypeListList{i}; % list of property types for this element
    propArrayList = ply.propArrayListList{i}; % list of property arrays for this element
    propCount = length(propNameList); % count of properties for this element
    dataCount = size(propArrayList{1},1); % length of data for this element
    
    % Write all data for element,
    % but process data differently for ascii and binary files.
    if strcmp(format,'ascii')
        
        % Put all property arrays for this element into a double array.
        A = zeros(dataCount,propCount);
        for j = 1:propCount
            A(:,j) = double(propArrayList{j});
        end
        
        % Construct format spec for fprintf.
        formatSpec = '';
        for j = 1:propCount
            switch propTypeList(j)
                case {'char', 'short', 'int'}
                    formatSpec = [formatSpec '%d'];
                case {'uchar', 'ushort', 'uint'}
                    formatSpec = [formatSpec '%u'];
                case {'float', 'double'}
                    formatSpec = [formatSpec '%g'];
                otherwise % include 'list'
                    error('property type "%s" not supported\n',propTypeList{j});
            end
            if j < propCount
                formatSpec = [formatSpec ' '];
            else
                formatSpec = [formatSpec '\n'];
            end
        end
        
        % Write all data for element.
        fprintf(fid,formatSpec,A');
        
    else % binary
        
        % Compute stride (length in bytes of all properties).
        stride = 0;
        for j = 1:propCount
            switch propTypeList(j)
                case {'char', 'uchar'}
                    stride = stride + 1;
                case {'short', 'ushort'}
                    stride = stride + 2;
                case {'int', 'uint', 'float'}
                    stride = stride + 4;
                case 'double'
                    stride = stride + 8;
                otherwise % include 'list'
                    error('property type "%s" not supported\n',propTypeList{j});
            end
        end
        
        % Process each property in turn.
        propFilePos = elemFilePos; % start at start of element
        for j = 1:propCount
            
            % Write property for entire element.
            switch propTypeList(j)
                case {'char', 'uchar'}
                    propBytes = 1;
                case {'short', 'ushort'}
                    propBytes = 2;
                case {'int', 'uint', 'float'}
                    propBytes = 4;
                case 'double'
                    propBytes = 8;
                otherwise % include 'list'
                    error('property type "%s" not supported\n',propTypeList(j));
            end
            if fseek(fid,propFilePos-(stride-propBytes),'bof') ~= 0
                error('element too wide compared to header');
            end
            fwrite(fid,propArrayList{j},propTypeList(j),stride-propBytes);
            propFilePos = propFilePos + propBytes;
        end
        elemFilePos = elemFilePos + (elemCount*stride); % file position of next element        
    end    
end

fclose(fid);
return;