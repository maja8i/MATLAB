%Display voxelized normals of reconstructed voxels by using the normal(s)
%of the nearest voxel(s) in the original (input) point cloud.
ptcloud_file_orig = '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10_WithNormalsAndCentroids\longdress_1300_voxelized10.ply';
ptcloud_file_recon_rootdir = '\\pandora\builds\test\Data\Compression\PLY\Codec_Results\longdress_1300\voxelized10\BV\';
ptcloud_recon_name = 'longdress_1300_voxelized10_distorted03';
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
%Sort the input normals (nx, ny, nz) in the above Morton order
normals_orig = ptcloud_orig(I, 4:6);
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

%SpatialIndex offset of an octree cell and its 26 neighbours
neighborOffset = [-1 -1 -1; -1 -1 0; -1 -1 1; -1 0 -1; -1 0 0; -1 0 1; -1 1 -1; -1 1 0; -1 1 1;
 0 -1 -1;  0 -1 0;  0 -1 1;  0 0 -1;  0 0 0;  0 0 1;  0 1 -1;  0 1 0;  0 1 1;
 1 -1 -1;  1 -1 0;  1 -1 1;  1 0 -1;  1 0 0;  1 0 1;  1 1 -1;  1 1 0;  1 1 1];

%Initialize a matrix to store the reconstructed voxel normals
normals_recon = zeros(size(xyz_recon, 1), 3);

tic;
%For each occupied octree cell at level lvl in the RECONSTRUCTED point 
%cloud ...
all_occ_nodes = (1:myOT_recon.NodeCount(b + 1))';
for lvl = (b + 1):-1:3
    normals_recon_cnt = 0;
    for occ_cell = all_occ_nodes'
        %Get the SpatialIndex of the current occ_cell
        spatialIndexOfBlock = double(myOT_recon.SpatialIndex{lvl}(occ_cell, :));
        %Get the SpatialIndex of each of the occ_cell's 26 neighbours
        spatialIndicesOfNeighboringBlocks = neighborOffset + repmat(spatialIndexOfBlock, 27, 1);
        spatialIndex = spatialIndicesOfNeighboringBlocks(all(spatialIndicesOfNeighboringBlocks >= 0, 2), :);
        %Get pointers to the chosen neighbourhood blocks (including the 
        %current occ_cell) in the ORIGINAL point cloud
        nP = myOT_orig.nodePtr(spatialIndex, (lvl - 1));
        %Get rid of neighbours that are not occupied (i.e., only keep 
        %occupied cells in the 27-neighbourhood)
        nP = nP(nP > 0); 
        if isempty(nP)
            if lvl == b + 1
                %If no corresponding neighbourhood is found in the original
                %point cloud, at the voxel level, mark the corresponding 
                %reconstructed voxels' normals as [NaN NaN NaN] for now
                normals_recon(((normals_recon_cnt + 1):(normals_recon_cnt + myOT_recon.DescendantCount{lvl}(occ_cell))), :) = NaN;
                normals_recon_cnt = normals_recon_cnt + myOT_recon.DescendantCount{lvl}(occ_cell);
            end
            continue;
        end
        %Count the total number of occupied voxels belonging to the current
        %occ_cell and its 26 neighbours in the ORIGINAL point cloud
        voxelCount = sum(myOT_orig.DescendantCount{lvl}(nP));
        %Allocate space to store the x, y, z coordinates of all the 
        %occupied voxels belonging to the octree cells in the current 
        %27-neighbourhood
        orig_voxel_subset = zeros(voxelCount, 3);
        %Allocate space for the normals of all the occupied voxels in
        %orig_voxel_subset
        voxel_subset_normals = zeros(voxelCount, 3);
        %Initialize counter for the total number of voxels (out of 
        %voxelCount) processed so far
        vox_count = 0;
        %Get the x, y, z coordinates of all the ORIGINAL occupied voxels 
        %belonging to the octree cells at level "lvl" denoted by nP, and 
        %their corresponding normals
        for n = nP'  
            %Get the index of the first occupied voxel in the n-th occupied 
            %cell at level lvl
            begin = myOT_orig.FirstDescendantPtr{lvl}(n);   
            %Get the total number of occupied voxels in the n-th occupied 
            %cell at level lvl
            count = myOT_orig.DescendantCount{lvl}(n);
            %Get the x, y, z location values for all the occupied voxels 
            %(at level b + 1) corresponding to the n-th occupied cell at 
            %level lvl
            orig_voxel_subset(((vox_count + 1):(vox_count + count)), :) = xyz_orig(begin:(begin + count - 1), :);
            %Get the nx, ny, nz normal values for the current occupied 
            %voxels
            voxel_subset_normals(((vox_count + 1):(vox_count + count)), :) = normals_orig(begin:(begin + count - 1), :);
            vox_count = vox_count + count;
        end

        %Get all the RECONSTRUCTED voxels belonging to the current occ_cell 
        %in the reconstructed point cloud
        begin2 = myOT_recon.FirstDescendantPtr{lvl}(occ_cell);
        count2 = myOT_recon.DescendantCount{lvl}(occ_cell);
        if lvl == b + 1
            recon_voxel_subset = xyz_recon(begin2:(begin2 + count2 - 1), :);
        else
            %Only get the reconstructed voxels that currently have NaN
            %normal values
            current_normals_recon_row_inds = begin2:(begin2 + count2 - 1);
            current_nan_normals_rows = all(isnan(normals_recon(current_normals_recon_row_inds, :)), 2);
            recon_voxel_subset = xyz_recon(current_normals_recon_row_inds(current_nan_normals_rows), :);
        end
        %For each of the reconstructed voxels in recon_voxel_subset, find 
        %their nearest neighbour inside orig_voxel_subset. Need to 
        %replicate orig_voxel_subset and recon_voxel_subset, so we can do a 
        %direct matrix subtraction between them.
        recon_voxel_subset_rep = recon_voxel_subset(repmat(1:size(recon_voxel_subset, 1), size(orig_voxel_subset, 1), 1), :);    %Replicate each row "size(orig_voxel_subset, 1)" times
        orig_voxel_subset_rep = repmat(orig_voxel_subset, size(recon_voxel_subset, 1), 1);    %Replicate entire matrix "size(recon_voxel_subset, 1)" times
        diffs = recon_voxel_subset_rep - orig_voxel_subset_rep;
        dists = sum(diffs.^2, 2);
        dists_reshaped = reshape(dists, size(recon_voxel_subset, 1), size(orig_voxel_subset, 1));
        min_dists = min(dists_reshaped, [], 2);
        %For each reconstructed voxel within the current occ_cell, if more 
        %than one equidistant original voxel is found, average the 
        %equidistant original voxels' normals (nx, ny, nz components 
        %separately) and map these average normal values to the 
        %corresponding reconstructed voxel
        recon_voxel_subset_normals = zeros(size(recon_voxel_subset));
        for rec_vox = 1:size(recon_voxel_subset, 1)
            recon_voxel_subset_normals(rec_vox, :) = mean(voxel_subset_normals((dists_reshaped(rec_vox, :) == min_dists(rec_vox)), :), 1);
        end
        if lvl == b + 1
            normals_recon(((normals_recon_cnt + 1):(normals_recon_cnt + count2)), :) = recon_voxel_subset_normals;
            normals_recon_cnt = normals_recon_cnt + count2;
        else
            normals_recon(current_normals_recon_row_inds(current_nan_normals_rows), :) = recon_voxel_subset_normals;
        end
    end %End occ_cell    
    %Check if there are any [NaN NaN NaN] normal values inside
    %normals_recon
    nan_rows = find(any(isnan(normals_recon), 2) == 1);
    %If yes, look at the next octree level up, to try to find nearest 
    %neighbours of the corresponding voxels 
    if ~isempty(nan_rows)
        %For each voxel indexed in nan_rows, find the pointer to its
        %ancestor at level "lvl - 1"
        firstDescendantPtr = myOT_recon.FirstDescendantPtr{lvl - 1};
        tmp = zeros(size(xyz_recon, 1), 1);
        tmp(firstDescendantPtr) = 1;
        ancestorPtrs_all = cumsum(tmp); %Ancestor cells at level "lvl - 1" for each voxel
        ancestorPtrs = ancestorPtrs_all(nan_rows);  %Extract only the ancestor pointers for voxels indexed in nan_rows
        all_occ_nodes = unique(ancestorPtrs);   %As some cells may have the same ancestor
    else
        break;
    end
end %End lvl

       
end_time = toc;
disp('************************************************************');
disp(['Time taken to find normals for reconstructed voxels: ' num2str(end_time) ' seconds']);
disp('************************************************************');

%Plot the reconstructed point cloud
a = [xyz_recon normals_recon];
figure; 
scatter3(a(:, 1), a(:, 2), a(:, 3), 5, [a(:, 4), a(:, 5), a(:, 6)], 'filled');
axis equal; axis off;

%Write the reconstructed point cloud with normals to PLY file ...

%Make a copy of plyStruct_recon, as it will be modified
plyStruct2 = plyStruct_recon;
%Create a new cell array for the plyStruct2 property ARRAYS, which contains
%the reconstructed voxel x, y, z coordinates and the normal nx, ny, nz data
plyStruct2.propArrayListList = cell(1, 1);
plyStruct2.propArrayListList{1}{1} = xyz_recon(:, 1);   %Reconstructed voxel X coordinates
plyStruct2.propArrayListList{1}{2} = xyz_recon(:, 2);   %Reconstructed voxel Y coordinates
plyStruct2.propArrayListList{1}{3} = xyz_recon(:, 3);   %Reconstructed voxel Z coordinates
plyStruct2.propArrayListList{1}{4} = normals_recon(:, 1);   %nx normal value
plyStruct2.propArrayListList{1}{5} = normals_recon(:, 2);   %ny normal value
plyStruct2.propArrayListList{1}{6} = normals_recon(:, 3);   %nz normal value
%Create a new cell array for the plyStruct2 property TYPES, which contains 
%the data types of the reconstructed voxel x, y, z coordinates, and the 
%data types of the normals
plyStruct2.propTypeListList = cell(1, 1);
plyStruct2.propTypeListList{1}(1) = "float";
plyStruct2.propTypeListList{1}(2) = "float";
plyStruct2.propTypeListList{1}(3) = "float";
plyStruct2.propTypeListList{1}(4) = "float";
plyStruct2.propTypeListList{1}(5) = "float";
plyStruct2.propTypeListList{1}(6) = "float";
%Create a new cell array for the plyStruct2 property NAMES, which contains 
%the names for the reconstructed voxel x, y, z coordinates, and the names
%for the nx, ny, nz components of the normals data 
plyStruct2.propNameListList = cell(1, 1);
plyStruct2.propNameListList{1}(1) = "x";
plyStruct2.propNameListList{1}(2) = "y";
plyStruct2.propNameListList{1}(3) = "z";
plyStruct2.propNameListList{1}(4) = "nx";
plyStruct2.propNameListList{1}(5) = "ny";
plyStruct2.propNameListList{1}(6) = "nz";

%Save the point cloud in the same directory as ptcloud_file_recon, but with
%"_normals" at the end of the name
plyWrite(plyStruct2, [ptcloud_file_recon '_normals.ply'], format);
