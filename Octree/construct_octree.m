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

%-------------------------------------------------------------------------%

function [myOT, mortonCodes_sorted, xyz_sorted, normals_sorted] = construct_octree(ptcloud_file, b)

%Read in input PLY point cloud
[~, ptcloud, ~] = plyRead(ptcloud_file);

%Get Morton codes for all x, y, z coordinates in the input point cloud
mortonCodes = xyzToMorton(ptcloud(:, 1:3), b);   %b bits for each Morton code
disp('Morton codes computed for input x, y, z');

%Sort the Morton codes obtained above, in ascending order
[mortonCodes_sorted, I] = sort(mortonCodes);
disp('Morton codes sorted');

%Arrange the input x, y, z locations in the same order as the sorted Morton
%codes, and the input normals in the same order (if the input point cloud
%has normals)
xyz_sorted = ptcloud(I, 1:3);
disp('x, y, z locations sorted');
if size(ptcloud, 2) == 9
    normals_sorted = ptcloud(I, 4:6);
    disp('Normals sorted');
end

%Construct an octree based on the sorted Morton codes computed above
disp(' ');
disp('Constructing octree ...');
myOT = octreeClass(mortonCodes_sorted, b);    %b-level octree
disp('Octree constructed');
disp('------------------------------------------------------------');


