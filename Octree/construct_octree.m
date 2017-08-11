%-------------------------------------------------------------------------%

%Construct an octree based on sorted Morton codes from an input 3D point
%cloud.

%Requires directory "Phil": add this directory plus its sub-directories to
%the current MATLAB path.

%---- INPUTS ----

%ptcloud_file: Name of the input PLY file, including the full path to this
%              file if it is not in the current MATLAB directory, e.g., 
%              ptcloud_file =
%              '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10\boxer_voxelized10.ply'.

%b: Bit depth for Morton codes and octree. b also determines the number of
%   levels in the octree that will be generated (apart from the root
%   level). The total number of octree levels INCLUDING the root level will
%   therefore be: b + 1. IMPORTANT: If using a voxelized point cloud as 
%   input, b must be equal to the voxelization level used to produce that 
%   point cloud, e.g., b = 10 for voxelized10 point clouds, b = 11 for 
%   voxelized11 clouds, etc.

%ptcloud_name: Just the name of the original (unvoxelized) input point 
%              cloud, without the PLY file extension, e.g., 
%              ptcloud_name = boxer   

%---- OUTPUTS ----

%myOT: Octree class, with a number of different properties.

%mortonCodes_sorted: Morton codes representing the (x, y, z) locations in
%                    the input point cloud, sorted in ascending order.

%xyz_sorted: Input (x, y, z) triplets arranged in the same order as their
%            corresponding Morton codes.

%Optional outputs:

%normals_sorted: Input normals (also (x, y, z) triplets) arranged in the
%                same order as their corresponding (x, y, z) location
%                values.

%centroids_sorted: Input centroids (also (x, y, z) triplets) arranged in 
%                  the same order as their corresponding (x, y, z) location 
%                  values.

%original_points_per_voxel: Original (x, y, z) points that were voxelized 
%                           to each voxel, if the input point cloud is a
%                           voxelized version of a non-voxelized point
%                           cloud.

%original_normals_per_voxel: The normal vectors (x, y, z) for each of the
%                            original_points_per_voxel, above.

%-------------------------------------------------------------------------%

function [myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, varargout] = construct_octree(ptcloud_file, b, varargin)

disp(['b = ' num2str(b)]);
disp(' ');

if numel(varargin) == 1
    ptcloud_name = varargin{1};
end

%Read in input PLY point cloud
[~, ptcloud, ~] = plyRead(ptcloud_file);

%Get Morton codes for all x, y, z coordinates in the input point cloud
mortonCodes = xyzToMorton(ptcloud(:, 1:3), b);   %b bits for each Morton code
disp('Morton codes computed for input x, y, z');

%Sort the Morton codes obtained above, in ascending order
[mortonCodes_sorted, I] = sort(mortonCodes);
disp('Morton codes sorted');

%Arrange the input x, y, z locations in the same order as the sorted Morton
%codes, and the input normals and centroids in the same order (if the input
%point cloud file contains normals and/or centroids)
xyz_sorted = ptcloud(I, 1:3);
disp('x, y, z locations sorted');
if size(ptcloud, 2) >= 9
    normals_sorted = ptcloud(I, 4:6);   %Assuming normals are in columns 4:6
    disp('Normals sorted');
end
if size(ptcloud, 2) == 12
    centroids_sorted = ptcloud(I, 10:12);   %Assuming centroids are in columns 10:12
    disp('Centroids sorted');
    varargout{1} = centroids_sorted;
end

%Construct an octree based on the sorted Morton codes computed above
disp(' ');
disp('Constructing octree ...');
myOT = octreeClass(mortonCodes_sorted, b);    %b-level octree
disp('Octree constructed');

if numel(varargin) == 1
    %If the input point cloud in ptcloud_file was a voxelized point cloud,
    %read in the original, unvoxelized version of that point cloud, and 
    %figure out, for each voxel, which (x, y, z) points in the original 
    %file were mapped to this voxel
    if b ~= 0
        disp(' ');
        %Read in the original, unvoxelized point cloud
        ptcloud_file_orig = ['\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\OriginalPLY_WithNormals\' ptcloud_name '.ply'];
        disp('Reading in corresponding unvoxelized point cloud ...');
        [~, ptcloud_orig, ~] = plyRead(ptcloud_file_orig);

        %Extract the (x, y, z) values of the points in ptcloud_orig
        xyz_orig = single(ptcloud_orig(:, 1:3));
        %Extract the normal (x, y, z) values for each point in xyz_orig
        if size(ptcloud_orig, 2) >= 9
            xyz_normals = single(ptcloud_orig(:, 7:9));
        end

        disp('Computing integer (voxel) coordinates ...');
        %Get the minimum x, y, and z values in xyz_orig
        minv = min(xyz_orig);
        %Get the maximum x, y, and z values in xyz_orig
        maxv = max(xyz_orig);
        %Get the maximum difference between any of the (maxv - minv) x, y, 
        %or z values
        width = max(maxv - minv);
        %Compute the integer location values for each original x, y, z 
        %value
        frame_to_world_scale = width/(2^b - 1);
        xyz_rounded = double(round((xyz_orig - repmat(minv, size(xyz_orig, 1), 1))/frame_to_world_scale));

        %Get Morton codes for all x, y, z coordinates in xyz_rounded
        mortonCodes2 = xyzToMorton(xyz_rounded, b);   %b bits for each Morton code
        disp('Morton codes computed for integer x, y, z');
        %Sort the Morton codes obtained above, in ascending order
        [~, I2] = sort(mortonCodes2);
        disp('Morton codes sorted');
        %Arrange the rows of xyz_rounded in the same order as their sorted
        %Morton codes above
        xyz_rounded_sorted = xyz_rounded(I2, 1:3);
        disp('Integer x, y, z sorted');

        %All the points that have been quantized to the same coordinate 
        %values belong to one voxel
        [xyz_rounded_sorted_unique, ~, voxel_inds] = unique(xyz_rounded_sorted, 'rows', 'stable');
        %Arrange the voxel indices in voxel_inds in the same order as their
        %corresponding xyz_rounded_sorted values
        voxel_inds_sorted = voxel_inds(I2);

        %Get all the original (x, y, z) points that belong to each voxel as
        %indicated by voxel_inds_sorted, and, if the input point cloud has
        %normals, then also get the corresponding normal vectors for the 
        %points in original_points_per_voxel, sorted in the same order
        original_points_per_voxel = cell(size(xyz_rounded_sorted_unique, 1), 1);
        if size(ptcloud_orig, 2) >= 9
            original_normals_per_voxel = cell(size(xyz_rounded_sorted_unique, 1), 1);
        end
        disp('Extracting original (x, y, z) points and normals that were mapped to each voxel ...');
        for vox = 1:size(xyz_rounded_sorted_unique, 1)
            original_points_per_voxel{vox} = xyz_orig((voxel_inds_sorted == vox), :);
            original_normals_per_voxel{vox} = xyz_normals((voxel_inds_sorted == vox), :);
            disp(['Finished voxel ' num2str(vox) '/' num2str(size(xyz_rounded_sorted_unique, 1))]);
        end
    end %End check if b ~= 0
end %End check if numel(varargin) == 1
    
%     %Extract the transform parameters
%     disp('Extracting transform parameters ...');
%     frame_to_world_translation = str2num(plyStruct.comments{3}(36:end));
%     frame_to_world_scale = str2double(plyStruct.comments{2}(30:end));
%     
%     %Translate xyz_orig to the coordinate system of the voxelized point 
%     %cloud
%     xyz_orig_transformed = xyz_orig - frame_to_world_translation;
%     %Scale xyz_orig_transformed
%     xyz_orig_transformed = xyz_orig_transformed./frame_to_world_scale;
%     %Quantize (round to nearest integer) xyz_orig_transformed
%     xyz_rounded = round(xyz_orig_transformed);
%     disp(' ');
%     disp('Transformed and quantized original (unvoxelized) x, y, z');
%        
%     %Get Morton codes for all x, y, z coordinates in xyz_rounded
%     mortonCodes2 = xyzToMorton(xyz_rounded, b);   %b bits for each Morton code
%     disp('Morton codes computed for transformed and quantized x, y, z');
%     %Sort the Morton codes obtained above, in ascending order
%     [mortonCodes2_sorted, I2] = sort(mortonCodes2);
%     disp('Morton codes sorted');
%     %Arrange the rows of xyz_rounded in the same order as their sorted
%     %Morton codes above
%     xyz_rounded_sorted = xyz_rounded(I2, 1:3);
%     disp('Transformed and quantized x, y, z sorted');
%     %All the points that have been quantized to the same coordinate values
%     %belong to one voxel
%     [xyz_rounded_sorted_unique, ~, voxel_inds] = unique(xyz_rounded_sorted, 'rows', 'stable');
%     %Arrange the voxel indices in voxel_inds in the same order as their
%     %corresponding xyz_rounded_sorted values
%     voxel_inds_sorted = voxel_inds(I2);
%     %Get all the original (x, y, z) points that belong to each voxel
%     %indicated by voxel_inds
%     original_points_per_voxel = xyz_orig(voxel_inds_sorted, :);
%     disp('Obtained original x, y, z points that were mapped to each voxel');
   
% end
disp('------------------------------------------------------------');

