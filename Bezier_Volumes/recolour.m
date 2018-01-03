%Colour reconstructed voxels by using the colour of the nearest voxel in
%the original (input) point cloud.
function [ptcloud_recon, mortonCodes_recon_sorted] = recolour(ptcloud_orig, ptcloud_recon, b, myOT_orig)

%Get the (Morton-sorted) input x, y, z locations and their corresponding R,
%G, B values
xyz_orig = ptcloud_orig(:, 1:3);
colours_orig = ptcloud_orig(:, 7:9);

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

%Initialize a matrix to store the reconstructed voxel colours
colours_recon = zeros(size(xyz_recon, 1), 3, 'single');

tic;
%For each occupied octree cell at level lvl in the RECONSTRUCTED point 
%cloud ...
all_occ_nodes = (1:myOT_recon.NodeCount(b + 1))';
for lvl = (b + 1):-1:3
    colours_recon_cnt = 0;
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
                %reconstructed voxels' colours as [NaN NaN NaN] for now
                colours_recon(((colours_recon_cnt + 1):(colours_recon_cnt + myOT_recon.DescendantCount{lvl}(occ_cell))), :) = NaN;
                colours_recon_cnt = colours_recon_cnt + myOT_recon.DescendantCount{lvl}(occ_cell);
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
        %Allocate space for the colours of all the occupied voxels in
        %orig_voxel_subset
        voxel_subset_colours = zeros(voxelCount, 3, 'single');
        %Initialize counter for the total number of voxels (out of 
        %voxelCount) processed so far
        vox_count = 0;
        %Get the x, y, z coordinates of all the ORIGINAL occupied voxels 
        %belonging to the octree cells at level "lvl" denoted by nP, and 
        %their corresponding colours
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
            %Get the R, G, B colour values for the current occupied voxels
            voxel_subset_colours(((vox_count + 1):(vox_count + count)), :) = colours_orig(begin:(begin + count - 1), :);
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
            %colour values
            current_colours_recon_row_inds = begin2:(begin2 + count2 - 1);
            current_nan_cols_rows = all(isnan(colours_recon(current_colours_recon_row_inds, :)), 2);
            recon_voxel_subset = xyz_recon(current_colours_recon_row_inds(current_nan_cols_rows), :);
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
        %equidistant original voxels' colours and map this average colour 
        %to the corresponding reconstructed voxel
        recon_voxel_subset_colours = zeros(size(recon_voxel_subset));
        for rec_vox = 1:size(recon_voxel_subset, 1)
            recon_voxel_subset_colours(rec_vox, :) = mean(voxel_subset_colours((dists_reshaped(rec_vox, :) == min_dists(rec_vox)), :), 1);
        end
        if lvl == b + 1
            colours_recon(((colours_recon_cnt + 1):(colours_recon_cnt + count2)), :) = uint8(recon_voxel_subset_colours);
            colours_recon_cnt = colours_recon_cnt + count2;
        else
            colours_recon(current_colours_recon_row_inds(current_nan_cols_rows), :) = uint8(recon_voxel_subset_colours);
        end
    end %End occ_cell    
    %Check if there are any [NaN NaN NaN] colour values inside
    %colours_recon
    nan_rows = find(any(isnan(colours_recon), 2) == 1);
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

%Consolidate the reconstructed x, y, z positions and their associated
%colours
ptcloud_recon = [xyz_recon colours_recon];
 
recolouring_time = toc;
disp('************************************************************');
disp(['Time taken to recolour reconstructed voxels: ' num2str(recolouring_time) ' seconds']);
disp('************************************************************');

% %Plot the reconstructed point cloud
% figure; 
% scatter3(ptcloud_recon(:, 1), ptcloud_recon(:, 2), ptcloud_recon(:, 3), 5, [ptcloud_recon(:, 4)./255, ptcloud_recon(:, 5)./255, ptcloud_recon(:, 6)./255], 'filled');
% axis equal; axis off;

