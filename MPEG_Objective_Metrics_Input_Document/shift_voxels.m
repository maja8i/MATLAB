%------------------------------ User Inputs ------------------------------%

%Below assumes that the input point cloud PLY file is in the current 
%directory
ptcloud_name = 'longdress_vox10_1330';    %Don't include ".ply" extension    
shifts = [1 2 3 4];
shift_dir = 0;  %0 => x; 1 => y; 2 =>z

%-------------------------------------------------------------------------% 

ptcloud_file = [ptcloud_name '.ply'];

%Read in input point cloud
[plyStruct, ptcloud, format] = plyRead(ptcloud_file);

%Extract just the voxel x, y, z coordinates
xyz = ptcloud(:, 1:3);

%Shift the input voxel x positions by each of the values in "shifts" in
%turn
for s = shifts
    %Shift all voxels by s in shift_dir direction
    if shift_dir == 0
        %Shift in x direction
        xyz2 = [(xyz(:, 1) + s) xyz(:, 2) xyz(:, 3)];
    elseif shift_dir == 1
        %Shift in y direction
        xyz2 = [xyz(:, 1) (xyz(:, 2) + s) xyz(:, 3)];
    elseif shift_dir == 2
        %Shift in z direction
        xyz2 = [xyz(:, 1) xyz(:, 2) (xyz(:, 3) + s)];
    end
    
    %Save the shifted point cloud as a PLY file (in the current directory)
    plyStruct2 = plyStruct;
    plyStruct2.propArrayListList = cell(1, 1);
    plyStruct2.propArrayListList{1}{1} = xyz2(:, 1);    %x   
    plyStruct2.propArrayListList{1}{2} = xyz2(:, 2);    %y  
    plyStruct2.propArrayListList{1}{3} = xyz2(:, 3);    %z   
    plyStruct2.propArrayListList{1}{4} = ptcloud(:, 4);   %R 
    plyStruct2.propArrayListList{1}{5} = ptcloud(:, 5);   %G
    plyStruct2.propArrayListList{1}{6} = ptcloud(:, 6);   %B
    
    plyStruct2.propTypeListList = cell(1, 1);
    plyStruct2.propTypeListList{1}(1) = "float";
    plyStruct2.propTypeListList{1}(2) = "float";
    plyStruct2.propTypeListList{1}(3) = "float";
    plyStruct2.propTypeListList{1}(4) = "uchar";
    plyStruct2.propTypeListList{1}(5) = "uchar";
    plyStruct2.propTypeListList{1}(6) = "uchar";
    
    plyStruct2.propNameListList = cell(1, 1);
    plyStruct2.propNameListList{1}(1) = "x";
    plyStruct2.propNameListList{1}(2) = "y";
    plyStruct2.propNameListList{1}(3) = "z";
    plyStruct2.propNameListList{1}(4) = "red";
    plyStruct2.propNameListList{1}(5) = "green";
    plyStruct2.propNameListList{1}(6) = "blue";
    
    if shift_dir == 0
        plyWrite(plyStruct2, [ptcloud_name '_Xshift' num2str(s) '.ply'], format);
    elseif shift_dir == 1
        plyWrite(plyStruct2, [ptcloud_name '_Yshift' num2str(s) '.ply'], format);
    elseif shift_dir == 2
        plyWrite(plyStruct2, [ptcloud_name '_Zshift' num2str(s) '.ply'], format);
    end
end
