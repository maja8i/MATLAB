%Colour reconstructed voxels by using the colour of the nearest voxel in
%the original (input) point cloud.
ptcloud_file_orig = '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10_WithNormalsAndCentroids\redandblack_1550_voxelized10.ply';
ptcloud_file_recon_rootdir = '\\pandora\builds\test\Data\Compression\PLY\Codec_Results\redandblack_1550\voxelized10\BezierVolume_thresh1\';
ptcloud_recon_name = 'redandblack_1550_voxelized10_distorted05';
ptcloud_file_recon = [ptcloud_file_recon_rootdir ptcloud_recon_name];
b = 10;

%Read in original and reconstructed point clouds
[~, ptcloud_orig, format] = plyRead(ptcloud_file_orig);
[plyStruct_recon, ptcloud_recon, ~] = plyRead([ptcloud_file_recon '.ply']);

%Get Morton codes for all x, y, z coordinates in the input point cloud
mortonCodes = xyzToMorton(ptcloud_orig(:, 1:3), b);   %b bits for each Morton code
%Sort the Morton codes obtained above, in ascending order
[mortonCodes_sorted, I] = sort(mortonCodes);
%Arrange the input x, y, z locations in the same order as the sorted Morton
%codes
xyz_orig = ptcloud_orig(I, 1:3);
%Sort the input colours in the above Morton order
colours_orig = ptcloud_orig(I, 7:9);
%Construct an octree based on the sorted Morton codes
myOT_orig = octreeClass(mortonCodes_sorted, b);    %b-level octree

%Get Morton codes for all x, y, z coordinates in the reconstructed point
%cloud
mortonCodes_recon = xyzToMorton(ptcloud_recon(:, 1:3), b);   %b bits for each Morton code
%Sort the Morton codes obtained above, in ascending order
[mortonCodes_recon_sorted, I2] = sort(mortonCodes_recon);
%Arrange the reconstructed x, y, z locations in the same order as the 
%sorted Morton codes
xyz_recon = ptcloud_recon(I2, 1:3);
%Construct an octree based on the sorted Morton codes
myOT_recon = octreeClass(mortonCodes_recon_sorted, b);    %b-level octree

%SpatialIndex offset of a voxel and its 26 neighbours
neighborOffset = [-1 -1 -1; -1 -1 0; -1 -1 1; -1 0 -1; -1 0 0; -1 0 1; -1 1 -1; -1 1 0; -1 1 1;
 0 -1 -1;  0 -1 0;  0 -1 1;  0 0 -1;  0 0 0;  0 0 1;  0 1 -1;  0 1 0;  0 1 1;
 1 -1 -1;  1 -1 0;  1 -1 1;  1 0 -1;  1 0 0;  1 0 1;  1 1 -1;  1 1 0;  1 1 1];

%Initialize a matrix to store the reconstructed voxel colours
colours_recon = zeros(size(ptcloud_recon, 1), 3, 'single');

tic;
%For each reconstructed voxel ...
for v = 1:size(ptcloud_recon(:, 1:3), 1)
    %Get the SpatialIndex of the current voxel
    spatialIndexOfBlock = double(myOT_recon.SpatialIndex{b + 1}(v, :));
    %Get the SpatialIndex of each of the current voxel's 26 neighbours
    spatialIndicesOfNeighboringBlocks = neighborOffset + repmat(spatialIndexOfBlock, 27, 1);
    spatialIndex = spatialIndicesOfNeighboringBlocks(all(spatialIndicesOfNeighboringBlocks >= 0, 2), :);
    %Get pointers to the chosen neighbourhood voxels (including the current
    %reconstructed voxel) in the original point cloud (if they exist)
    nP = myOT_orig.nodePtr(spatialIndex, b);
    %Get rid of neighbours that aren't occupied (i.e., only keep occupied
    %neighbours)
    nP = nP(nP > 0); 
    %Get the x, y, z coordinates of the original voxels indexed by nP
    orig_voxel_subset = xyz_orig(nP, :);
    %Get the colours (R, G, B values) of the original voxels indexed by nP
    orig_voxel_colours = colours_orig(nP, :);
    %Find the nearest voxel in orig_voxel_subset, to the current
    %reconstructed voxel
    diffs = ptcloud_recon(v, :) - orig_voxel_subset;
    dists = sum(diffs.^2, 2);
    min_dist = min(dists);
    %If more than one equidistant voxel is found, average the equidistant
    %voxels' colours and map this average colour to the current
    %reconstructed voxel
    colours_recon(v, :) = mean(orig_voxel_colours((dists == min_dist), :), 1);
end
recolouring_time = toc;
disp('************************************************************');
disp(['Time taken to recolour reconstructed voxels: ' num2str(recolouring_time) ' seconds']);
disp('************************************************************');

%Plot the reconstructed point cloud
a = [xyz_recon colours_recon];
figure; 
scatter3(a(:, 1), a(:, 2), a(:, 3), 5, [a(:, 4)./255, a(:, 5)./255, a(:, 6)./255], 'filled');
axis equal; axis off;

%Write the reconstructed point cloud with colours to PLY file ...

%Make a copy of plyStruct_recon, as it will be modified
plyStruct2 = plyStruct_recon;
%Create a new cell array for the plyStruct2 property ARRAYS, which contains
%the reconstructed voxel x, y, z coordinates and the colour R, G, B data
plyStruct2.propArrayListList = cell(1, 1);
plyStruct2.propArrayListList{1}{1} = xyz_recon(:, 1);   %Reconstructed voxel X coordinates
plyStruct2.propArrayListList{1}{2} = xyz_recon(:, 2);   %Reconstructed voxel Y coordinates
plyStruct2.propArrayListList{1}{3} = xyz_recon(:, 3);   %Reconstructed voxel Z coordinates
plyStruct2.propArrayListList{1}{4} = uint8(colours_recon(:, 1));   %R colour value
plyStruct2.propArrayListList{1}{5} = uint8(colours_recon(:, 2));   %G colour value
plyStruct2.propArrayListList{1}{6} = uint8(colours_recon(:, 3));   %B colour value
%Create a new cell array for the plyStruct2 property TYPES, which contains 
%the data types of the reconstructed voxel x, y, z coordinates, and the 
%data types of the colours 
plyStruct2.propTypeListList = cell(1, 1);
plyStruct2.propTypeListList{1}(1) = "float";
plyStruct2.propTypeListList{1}(2) = "float";
plyStruct2.propTypeListList{1}(3) = "float";
plyStruct2.propTypeListList{1}(4) = "uchar";
plyStruct2.propTypeListList{1}(5) = "uchar";
plyStruct2.propTypeListList{1}(6) = "uchar";
%Create a new cell array for the plyStruct2 property NAMES, which contains 
%the names for the reconstructed voxel x, y, z coordinates, and the names
%for the R, G, B channels of the colour data 
plyStruct2.propNameListList = cell(1, 1);
plyStruct2.propNameListList{1}(1) = "x";
plyStruct2.propNameListList{1}(2) = "y";
plyStruct2.propNameListList{1}(3) = "z";
plyStruct2.propNameListList{1}(4) = "red";
plyStruct2.propNameListList{1}(5) = "green";
plyStruct2.propNameListList{1}(6) = "blue";

%Save the point cloud in the same directory as ptcloud_file_recon, but with
%"_coloured" at the end of the name
plyWrite(plyStruct2, [ptcloud_file_recon '_coloured.ply'], format);
