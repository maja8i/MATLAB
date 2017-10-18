% Read generic PLY file
% (arbitrary element and property names, arbitrary types; list types not implemented)
%function plyStruct = plyRead(filename)
function [plyStruct, A, format] = plyRead(filename)     %Modified by Maja
%
% Input:
% filename - character string
%
% Output:
% plyStruct.elemNameList - elemCount x 1 string array of element names,
%   e.g., ["vertex"; "face"], where "vertex" means string('vertex')
% plyStruct.propNameListList - elemCount x 1 cell array of string arrays of property names,
%   e.g., {["x","y","z","red","green","blue"];["vertex_indices"]}
% plyStruct.propTypeListList - elemCount x 1 cell array of string arrays of property types,
%   e.g., {["float","float","float","uchar","uchar","uchar"};["list uchar int"]}
% plyStruct.propArrayListList - elemCount x 1 cell array of cell arrays of property arrays,
%   e.g., {{x,y,z,red,green,blue};{vertex_indices}}
%   where x, y, and z are nx1 float; red, green, and blue are nx1 uchar;
%   and vertex_indices are mxf int, where n and m are vertex and face counts
% plyStruct.comments - k x 1 string array of comment lines
%
% "list" types are currently unsupported

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

% Get the header into a cell array of character strings.
% Separate out the comments into their own cell array of character strings.
header = strings(0);
comments = strings(0);
fid = fopen(filename,'r','n','US-ASCII');
if fid < 0
    error('cannot open file "%s" for reading\n',filename);
end
s = string(fgetl(fid));
while strcmp(s,'end_header') == 0
    if startsWith(s,'comment')
        comments = [comments; s];
    else
        header = [header; s];
    end
    s = string(fgetl(fid));
end
headerByteCount = ftell(fid);
headerLineCount = length(header);
fclose(fid);

% Make sure it's a ply file.
if headerLineCount < 3 || strcmp(header(1),'ply') == 0
    error('incomplete ply file');
end

% Is it ascii or binary?
format = sscanf(header(2),'%*s %s %*s');

% Re-open file with the correct format.
switch format
    case 'ascii'
        fid = fopen(filename,'r','n','US-ASCII');
    case 'binary_little_endian'
        fid = fopen(filename,'r','ieee-le');
    case 'binary_big_endian'
        fid = fopen(filename,'r','ieee-be');
    otherwise
        error('format "%s" not supported\n',format);
end

% Seek to beginning of data.
if fseek(fid,headerByteCount,'bof') ~= 0
    error('incomplete ply file');
end

% Initialize lists to return.
elemNameList = strings(0); % List of element names
propArrayListList = cell(0); % List of lists of array containing properties
propNameListList = cell(0); % List of lists of associated property names
propTypeListList = cell(0); % List of lists of associated property types

% Process element-by-element.
elemFilePos = headerByteCount; % position of current element's data in file
elemLineNum = 3; % first element should be described on header line 3
while elemLineNum <= headerLineCount && startsWith(header(elemLineNum),'element')
    
    % Get element name and count.
    elemName = sscanf(header(elemLineNum),'%*s %s');
    elemCount = sscanf(header(elemLineNum),'%*s %*s %d');
    
    % Store element name (ignoring those with count==0).
    if elemCount > 0
        elemNameList = [elemNameList; elemName];
    end
    
    % Allocate array of length 'count' for each property.
    propNameList = strings(0);
    propTypeList = strings(0);
    propLineNum = elemLineNum + 1;
    while propLineNum <= headerLineCount && startsWith(header(propLineNum),'property')
        
        % Get property type and name.
        propSplit = split(header(propLineNum));
        propType = propSplit(2); % float, uchar, list, etc.
        if ~strcmp(propType,'list') % float, uchar, etc.
            propName = propSplit(3); % x, y, z, red, green, blue, etc.
        else % list
            listCountType = propSplit(3); % uchar, int, etc.
            listType = propSplit(4); % int, etc.
            propName = propSplit(5); % vertex_indices, etc.
        end
        
        % Save property type and name.
        propTypeList = [propTypeList, propType];
        propNameList = [propNameList, propName];
        
        % Get next property.
        propLineNum = propLineNum + 1;
    end
    
    propCount = propLineNum - elemLineNum - 1; % number of properties for this element
    propArrayList = cell(1,propCount); % cell array to hold array for each property
    
    % Don't read elements with count==0.
    if elemCount == 0
        elemLineNum = propLineNum;
        continue;
    end
    
    % Read all data for element,
    % but process data differently for ascii and binary files.
    if strcmp(format,'ascii')
        
        % Construct format spec for fscanf.
        formatSpec = '';
        for j = 1:propCount
            switch propTypeList(j)
                case {'char', 'short', 'int'}
                    formatSpec = [formatSpec '%d'];
                case {'uchar', 'ushort', 'uint'}
                    formatSpec = [formatSpec '%u'];
                case {'float', 'double'}
                    formatSpec = [formatSpec '%f'];
                otherwise % include 'list'
                    error('property type "%s" not supported\n',propTypeList{j});
            end
            if j < propCount
                formatSpec = [formatSpec ' '];
            else
                formatSpec = [formatSpec '\n'];
            end
        end
        
        % Read all data for element.
        A = fscanf(fid,formatSpec,[propCount,elemCount])';

        % Move it to arrays of the right type.
        for j = 1:propCount
            switch propTypeList(j)
                case 'char'
                    propArrayList{j} = int8(A(:,j));
                case 'uchar'
                    propArrayList{j} = uint8(A(:,j));
                case 'short'
                    propArrayList{j} = int16(A(:,j));
                case 'ushort'
                    propArrayList{j} = uint16(A(:,j));
                case 'int'
                    propArrayList{j} = int32(A(:,j));
                case 'uint'
                    propArrayList{j} = uint32(A(:,j));
                case 'float'
                    propArrayList{j} = single(A(:,j));
                case 'double'
                    propArrayList{j} = double(A(:,j));
            end
        end
        
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
            
            % Seek to beginning of property.
            if fseek(fid,propFilePos,'bof') ~= 0
                error('incomplete ply file');
            end
            
            % Read property for entire element.
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
            propArrayList{j} = fread(fid,[elemCount,1],['*' propTypeList{j}],stride-propBytes);
            propFilePos = propFilePos + propBytes;
        end
        elemFilePos = elemFilePos + (elemCount*stride); % file position of next element
    end
        
    % Save property lists.
    propArrayListList = cat(1, propArrayListList, {propArrayList});
    propNameListList = cat(1, propNameListList, {propNameList});
    propTypeListList = cat(1, propTypeListList, {propTypeList});
    
    % Get next element.
    elemLineNum = propLineNum;
end

% Should have processed all lines in header.
if elemLineNum <= headerLineCount
    disp('invalid ply header');
    return;
end

% Close file.
fclose(fid);

plyStruct = struct('elemNameList',elemNameList,...
    'propNameListList',{propNameListList},...
    'propTypeListList',{propTypeListList},...
    'propArrayListList',{propArrayListList},...
    'comments',comments);