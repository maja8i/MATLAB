function [reconstructed_vox_pos, reconstructed_vox_pos_corners] = voxel_reconstruction_pruning(pp_first_nonempty, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, ptcloud_file, ptcloud_name)

tic;

%Matrix to store the reconstructed voxel corners: initialize to some large
%value (can shrink later)
reconstructed_vox_pos_corners = zeros(24000000, 3);
%Counter for reconstructed_vox_pos_corners
rec_vox_pos_cnrs_cntr = 1;

%If there are any voxels that have already been reconstructed (because no
%pruning was done further up the octree for these voxels), get their corner
%coordinates
if ~isempty(corner_coords_decoder{b + 1})
    disp(' ');
    disp('Getting corner coordinates for already-reconstructed occupied voxels ...');
    disp(' ');    
    reconstructed_vox_pos_corners(rec_vox_pos_cnrs_cntr:size(corner_coords_decoder{b + 1}, 1), 1:3) = corner_coords_decoder{b + 1};
    rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + size(corner_coords_decoder{b + 1}, 1);
end

%Cell array that will store the indices of only the occupied cells from
%post_pruning_array that are leaves (the leaf cells are the only ones that
%may be subdivided further and serve as starting points for the trilinear 
%interpolation, not the internal cells, so we do not need to consider the
%internal cells for voxel reconstruction)
pp_leaves_only = cell(b, 1);
for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    %Find the location(s) of the leaf cell(s)
    leaves = find(post_pruning_array{lvl} == 1);
    %Store the indices of the leaf cells at this level
    pp_leaves_only{lvl} = leaves;
end

%For debugging purposes, print out how many (if any) of the leaf cells at 
%each octree level in pp_leaves_only, end up having all 8 of their 
%decoder-reconstructed control points with the same sign (+/-) ...

for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    same_sign_cntr = 0;
    zero_cp_cntr = 0;
    for occ_cell = pp_leaves_only{lvl}'
        %Get all 8 control points for the corners of the current leaf cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Check if all control points of the current leaf cell have the same
        %sign, including the case where all the control points may be 0
        if (abs(sum(sign(current_ctrlpts))) == 8)||(~any(sign(current_ctrlpts)))
            same_sign_cntr = same_sign_cntr + 1;
            if ~any(sign(current_ctrlpts))
                zero_cp_cntr = zero_cp_cntr + 1;
            end
            %disp(['Leaf cell ' num2str(occ_cell) ' has all control points with the same sign: ']);
            %Display the control points for all corners of this leaf cell
            %disp(num2str(current_ctrlpts));
            %disp(' ');
            %Get the corner coordinates for each corner of this leaf cell
            %cell_corners = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
            %disp('Leaf cell corner coordinates:');
            %isp(num2str(cell_corners));
            %disp(' ');
        end
    end %End i
    disp(['TOTAL number of leaf cells with all control points having the same sign, at level ' num2str(lvl) ': ' num2str(same_sign_cntr) '/' num2str(length(pp_leaves_only{lvl})) ' (' num2str((same_sign_cntr/(length(pp_leaves_only{lvl}))*100)) '%)']);
    disp(['No. of leaf cells with all 0 control points at level ' num2str(lvl) ': ' num2str(zero_cp_cntr)]);
    disp(' ');
end %End lvl

%For each octree level at which there are leaf cells, except the voxel
%level ...
for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    %profile on
    %For each leaf cell at this level ...
    for occ_cell = pp_leaves_only{lvl}'
        disp(['Reconstructing voxels for leaf cell (occ_cell) ' num2str(occ_cell) ' at level ' num2str(lvl) ':']);
        %Get the 8 corner coordinates of the current leaf cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
        %Get the control points for all 8 corners of the current leaf cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Keep subdividing the current leaf cell until we reach cells of 
        %size 1 x 1 x 1 (i.e., voxels). For each sub-cell at each level,
        %interpolate between the leaf cell's (occ_cell's) control points.
        %The sub-cells that have interpolated control points with different 
        %signs will be subdivided further; at the end, voxels that have 
        %interpolated control points with different signs will be 
        %considered occupied.
        for lvl_d = (lvl + 1):(b + 1)
            %Compute the corner coordinates of each of the sub-cells at
            %lvl_d ...
            if lvl_d == (lvl + 1)
                disp(['--Processing sub-cells at lvl_d ' num2str(lvl_d) ': 8 sub-cells to process at current lvl_d']);
                %If we are at the first level after the leaf, get the 
                %minimum, maximum and midpoint of the x, y, and z 
                %coordinates of the leaf cell. These will be used to 
                %compute the corner coordinates of the 8 sub-cells at level
                %lvl_d.  
                min_x = min(current_corner_coords(:, 1));  
                min_y = min(current_corner_coords(:, 2));
                min_z = min(current_corner_coords(:, 3));
                max_x = max(current_corner_coords(:, 1));
                max_y = max(current_corner_coords(:, 2));
                max_z = max(current_corner_coords(:, 3));
                mid_x = (min_x + max_x)/2;
                mid_y = (min_y + max_y)/2;
                mid_z = (min_z + max_z)/2;   
                zc_coords = current_corner_coords;  %"zc" stands for zero crossing, as we will only subdivide occupied cells
            else
                disp(['--Processing sub-cells at lvl_d ' num2str(lvl_d) ': ' num2str(size(subcell_coords_occupied, 1)/8) ' occupied sub-cells at previous lvl_d, so ' num2str(size(subcell_coords_occupied, 1)) ' sub-cells to process at current lvl_d']);    
                %If none of the sub-cells at the previous lvl_d were marked
                %as candidates for further subdivision (i.e., none of them
                %were considered occupied), then do not check any more sub-
                %cells for the current leaf cell (occ_cell), but move on to
                %checking the next leaf cell
                if isempty(subcell_coords_occupied)
                    break;
                end
                %Find the minimum, maximum and midpoint x, y, z of the 
                %corner coordinates of each of the cells at the previous 
                %level that have been marked as candidates for further 
                %subdivision (those stored in subcell_coords_occupied). 
                %There will be one min, one max and one mid point (for x, 
                %y, z separately) per 8 sub-cells at the current lvl_d 
                %(there are 8 sub-cells, not all necessarily occupied, per 
                %occupied cell in subcell_coords_occupied).
                min_x(1:(size(subcell_coords_occupied, 1)/8), 1) = min(reshape(subcell_coords_occupied(:, 1), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                min_y(1:(size(subcell_coords_occupied, 1)/8), 1) = min(reshape(subcell_coords_occupied(:, 2), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                min_z(1:(size(subcell_coords_occupied, 1)/8), 1) = min(reshape(subcell_coords_occupied(:, 3), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                max_x(1:(size(subcell_coords_occupied, 1)/8), 1) = max(reshape(subcell_coords_occupied(:, 1), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                max_y(1:(size(subcell_coords_occupied, 1)/8), 1) = max(reshape(subcell_coords_occupied(:, 2), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                max_z(1:(size(subcell_coords_occupied, 1)/8), 1) = max(reshape(subcell_coords_occupied(:, 3), 8, size(subcell_coords_occupied, 1)/8), [], 1);
                mid_x(1:length(min_x), 1) = (min_x + max_x)./2;
                mid_y(1:length(min_x), 1) = (min_y + max_y)./2;
                mid_z(1:length(min_x), 1) = (min_z + max_z)./2;
                zc_coords = subcell_coords_occupied;
            end %End check if lvl_d == (lvl + 1) 
            %Use the min, max, and mid points found above, to compute the
            %corner coordinates of each of the 8 sub-cells resulting from
            %subdividing each of the cells whose coordinates are stored in
            %zc_coords
            disp(['  Computing corner coordinates for each sub-cell at this level (' num2str(size(zc_coords, 1)) ' sub-cells in total)']);           
            %Corner coordinates of sub-cells 1 of each cell in zc_coords: 
            %each set of length(min_x) rows in subcell_coords_temp, below, 
            %represents one of the 8 corners of sub-cell 1 of each cell in 
            %zc_coords
            subcell_coords_temp = [min_x min_y min_z; mid_x min_y min_z; mid_x mid_y min_z; min_x mid_y min_z; min_x min_y mid_z; mid_x min_y mid_z; mid_x mid_y mid_z; min_x mid_y mid_z]; 
            %Rearrange subcell_coords_temp so that all 8 corners of one
            %sub-cell are one after the other
            cnrs1 = subcell_coords_temp(1:length(min_x), :);
            cnrs2 = subcell_coords_temp((length(min_x) + 1):(2*length(min_x)), :);
            cnrs3 = subcell_coords_temp((2*length(min_x) + 1):(3*length(min_x)), :);
            cnrs4 = subcell_coords_temp((3*length(min_x) + 1):(4*length(min_x)), :);
            cnrs5 = subcell_coords_temp((4*length(min_x) + 1):(5*length(min_x)), :);
            cnrs6 = subcell_coords_temp((5*length(min_x) + 1):(6*length(min_x)), :);
            cnrs7 = subcell_coords_temp((6*length(min_x) + 1):(7*length(min_x)), :);
            cnrs8 = subcell_coords_temp((7*length(min_x) + 1):(8*length(min_x)), :);
            subcell_coords_temp = zeros(size(subcell_coords_temp, 1), 3);
            subcell_coords_temp(1:8:(size(subcell_coords_temp, 1) - 7), 1:3) = cnrs1; 
            subcell_coords_temp(2:8:(size(subcell_coords_temp, 1) - 6), 1:3) = cnrs2; 
            subcell_coords_temp(3:8:(size(subcell_coords_temp, 1) - 5), 1:3) = cnrs3; 
            subcell_coords_temp(4:8:(size(subcell_coords_temp, 1) - 4), 1:3) = cnrs4; 
            subcell_coords_temp(5:8:(size(subcell_coords_temp, 1) - 3), 1:3) = cnrs5; 
            subcell_coords_temp(6:8:(size(subcell_coords_temp, 1) - 2), 1:3) = cnrs6; 
            subcell_coords_temp(7:8:(size(subcell_coords_temp, 1) - 1), 1:3) = cnrs7; 
            subcell_coords_temp(8:8:size(subcell_coords_temp, 1), 1:3) = cnrs8;  
            %Offsets for computing corner coordinates of sub-cells 1-8 for 
            %each cell in zc_coords (sub-cells 1 have already been computed
            %above, so the offsets for them are just [0 0 0])
            offsets = [0 0 0;
                2^(b + 1 - lvl_d) 0 0;
                2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d) 0;
                0 2^(b + 1 - lvl_d) 0;
                0 0 2^(b + 1 - lvl_d);
                2^(b + 1 - lvl_d) 0 2^(b + 1 - lvl_d);
                2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d);
                0 2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d)]; 
            %Corner coordinates of sub-cells 1-8 for each cell in zc_coords
            subcell_coords = subcell_coords_temp(repmat(1:size(subcell_coords_temp, 1), 8, 1), :) + repmat(offsets, size(zc_coords, 1), 1);
            %Normalize all subcell_coords to be in the range [0, 1],
            %because the trilinear interpolation formula (used below, to
            %obtain the interpolated control points) will expect the x, y,
            %z values to be in this range 
            subcell_coords_orig = subcell_coords;
            subcell_coords = (subcell_coords - current_corner_coords(1, :))/(2^(b + 1 - lvl));    %Subtract the origin of the original parent (leaf) cell, and divide by the parent cell width
            %Array that will store the corner coordinates of sub-cells at
            %the current lvl_d, which are candidates for further
            %subdivision (if lvl_d < b + 1)
            subcell_coords_occupied = [];
            %For each corner computed in subcell_coords, calculate the
            %interpolated control point for this corner, by trilinearly
            %interpolating between the original parent control points,
            %current_ctrlpts ...
            disp(['  Computing interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
            %Collect all "corner 1"s of the sub-cells inside subcell_coords
            corners1 = subcell_coords(1:8:(end - 7), :);
            %Collect all "corner 2"s of the sub-cells inside subcell_coords
            corners2 = subcell_coords(2:8:(end - 6), :);
            %Collect all "corner 3"s of the sub-cells inside subcell_coords
            corners3 = subcell_coords(3:8:(end - 5), :);
            %Collect all "corner 4"s of the sub-cells inside subcell_coords
            corners4 = subcell_coords(4:8:(end - 4), :);
            %Collect all "corner 5"s of the sub-cells inside subcell_coords
            corners5 = subcell_coords(5:8:(end - 3), :);
            %Collect all "corner 6"s of the sub-cells inside subcell_coords
            corners6 = subcell_coords(6:8:(end - 2), :);
            %Collect all "corner 7"s of the sub-cells inside subcell_coords
            corners7 = subcell_coords(7:8:(end - 1), :);
            %Collect all "corner 8"s of the sub-cells inside subcell_coords
            corners8 = subcell_coords(8:8:end, :);
            %Compute the multiplication matrix of corner coordinates for
            %all the sub-cells represented in subcell_coords: these will
            %each be multiplied by the control point of the corresponding
            %leaf cell corner, when doing the trilinear interpolation
            mult_matrix = [(1 - corners1(:, 1)).*(1 - corners1(:, 2)).*(1 - corners1(:, 3));
                    corners1(:, 1).*(1 - corners1(:, 2)).*(1 - corners1(:, 3));
                    corners1(:, 1).*corners1(:, 2).*(1 - corners1(:, 3));
                    (1 - corners1(:, 1)).*corners1(:, 2).*(1 - corners1(:, 3));
                    (1 - corners1(:, 1)).*(1 - corners1(:, 2)).*corners1(:, 3);
                    corners1(:, 1).*(1 - corners1(:, 2)).*corners1(:, 3);
                    corners1(:, 1).*corners1(:, 2).*corners1(:, 3);
                    (1 - corners1(:, 1)).*corners1(:, 2).*corners1(:, 3);  %End corners 1
                    (1 - corners2(:, 1)).*(1 - corners2(:, 2)).*(1 - corners2(:, 3));
                    corners2(:, 1).*(1 - corners2(:, 2)).*(1 - corners2(:, 3));
                    corners2(:, 1).*corners2(:, 2).*(1 - corners2(:, 3));
                    (1 - corners2(:, 1)).*corners2(:, 2).*(1 - corners2(:, 3));
                    (1 - corners2(:, 1)).*(1 - corners2(:, 2)).*corners2(:, 3);
                    corners2(:, 1).*(1 - corners2(:, 2)).*corners2(:, 3);
                    corners2(:, 1).*corners2(:, 2).*corners2(:, 3);
                    (1 - corners2(:, 1)).*corners2(:, 2).*corners2(:, 3);  %End corners 2
                    (1 - corners3(:, 1)).*(1 - corners3(:, 2)).*(1 - corners3(:, 3));
                    corners3(:, 1).*(1 - corners3(:, 2)).*(1 - corners3(:, 3));
                    corners3(:, 1).*corners3(:, 2).*(1 - corners3(:, 3));
                    (1 - corners3(:, 1)).*corners3(:, 2).*(1 - corners3(:, 3));
                    (1 - corners3(:, 1)).*(1 - corners3(:, 2)).*corners3(:, 3);
                    corners3(:, 1).*(1 - corners3(:, 2)).*corners3(:, 3);
                    corners3(:, 1).*corners3(:, 2).*corners3(:, 3);
                    (1 - corners3(:, 1)).*corners3(:, 2).*corners3(:, 3);  %End corners 3
                    (1 - corners4(:, 1)).*(1 - corners4(:, 2)).*(1 - corners4(:, 3));
                    corners4(:, 1).*(1 - corners4(:, 2)).*(1 - corners4(:, 3));
                    corners4(:, 1).*corners4(:, 2).*(1 - corners4(:, 3));
                    (1 - corners4(:, 1)).*corners4(:, 2).*(1 - corners4(:, 3));
                    (1 - corners4(:, 1)).*(1 - corners4(:, 2)).*corners4(:, 3);
                    corners4(:, 1).*(1 - corners4(:, 2)).*corners4(:, 3);
                    corners4(:, 1).*corners4(:, 2).*corners4(:, 3);
                    (1 - corners4(:, 1)).*corners4(:, 2).*corners4(:, 3);  %End corners 4
                    (1 - corners5(:, 1)).*(1 - corners5(:, 2)).*(1 - corners5(:, 3));
                    corners5(:, 1).*(1 - corners5(:, 2)).*(1 - corners5(:, 3));
                    corners5(:, 1).*corners5(:, 2).*(1 - corners5(:, 3));
                    (1 - corners5(:, 1)).*corners5(:, 2).*(1 - corners5(:, 3));
                    (1 - corners5(:, 1)).*(1 - corners5(:, 2)).*corners5(:, 3);
                    corners5(:, 1).*(1 - corners5(:, 2)).*corners5(:, 3);
                    corners5(:, 1).*corners5(:, 2).*corners5(:, 3);
                    (1 - corners5(:, 1)).*corners5(:, 2).*corners5(:, 3);  %End corners 5
                    (1 - corners6(:, 1)).*(1 - corners6(:, 2)).*(1 - corners6(:, 3));
                    corners6(:, 1).*(1 - corners6(:, 2)).*(1 - corners6(:, 3));
                    corners6(:, 1).*corners6(:, 2).*(1 - corners6(:, 3));
                    (1 - corners6(:, 1)).*corners6(:, 2).*(1 - corners6(:, 3));
                    (1 - corners6(:, 1)).*(1 - corners6(:, 2)).*corners6(:, 3);
                    corners6(:, 1).*(1 - corners6(:, 2)).*corners6(:, 3);
                    corners6(:, 1).*corners6(:, 2).*corners6(:, 3);
                    (1 - corners6(:, 1)).*corners6(:, 2).*corners6(:, 3);  %End corners 6
                    (1 - corners7(:, 1)).*(1 - corners7(:, 2)).*(1 - corners7(:, 3));
                    corners7(:, 1).*(1 - corners7(:, 2)).*(1 - corners7(:, 3));
                    corners7(:, 1).*corners7(:, 2).*(1 - corners7(:, 3));
                    (1 - corners7(:, 1)).*corners7(:, 2).*(1 - corners7(:, 3));
                    (1 - corners7(:, 1)).*(1 - corners7(:, 2)).*corners7(:, 3);
                    corners7(:, 1).*(1 - corners7(:, 2)).*corners7(:, 3);
                    corners7(:, 1).*corners7(:, 2).*corners7(:, 3);
                    (1 - corners7(:, 1)).*corners7(:, 2).*corners7(:, 3);  %End corners 7
                    (1 - corners8(:, 1)).*(1 - corners8(:, 2)).*(1 - corners8(:, 3));
                    corners8(:, 1).*(1 - corners8(:, 2)).*(1 - corners8(:, 3));
                    corners8(:, 1).*corners8(:, 2).*(1 - corners8(:, 3));
                    (1 - corners8(:, 1)).*corners8(:, 2).*(1 - corners8(:, 3));
                    (1 - corners8(:, 1)).*(1 - corners8(:, 2)).*corners8(:, 3);
                    corners8(:, 1).*(1 - corners8(:, 2)).*corners8(:, 3);
                    corners8(:, 1).*corners8(:, 2).*corners8(:, 3);
                    (1 - corners8(:, 1)).*corners8(:, 2).*corners8(:, 3)]; %End corners 8
            %Multiply each control point of the current leaf cell (parent
            %cell of the current sub-cells) by the corresponding lines in
            %mult_matrix
            current_ctrlpts_rep = current_ctrlpts(repmat(1:length(current_ctrlpts), (size(subcell_coords, 1)/8), 1), :);
            temp = repmat(current_ctrlpts_rep, 8, 1).*mult_matrix;
            %Compute the trilinear interpolation sum for each of the 8 
            %corners of each sub-cell
            temp_reshaped1 = reshape(temp, size(subcell_coords, 1), 8);
            subcell_ctrlpts_temp = zeros(size(subcell_coords, 1), 1);
            for cnr = 1:8
                %Reshape column "cnr" of temp_reshaped1 (this column
                %represents all corners "cnr" for all current sub-cells) so
                %that each column represents one multiplication from "temp"
                %(there are 8 multiplications per corner, so 8 columns) and
                %each row represents one sub-cell. Then sum up the
                %multiplications in each row of temp_reshaped2, and the
                %result will represent the interpolated control point for 
                %corner "cnr" of the sub-cell represented by that row of
                %temp_reshaped2.
                temp_reshaped2 = reshape(temp_reshaped1(:, cnr), (size(subcell_coords, 1)/8), 8);   %Each row represents one sub-cell, each column one of the 8 multiplications for corner "cnr"
                %Compute interpolated control point for corner "cnr" of 
                %each sub-cell
                subcell_ctrlpts_temp(((cnr*(size(subcell_coords, 1)/8) - (size(subcell_coords, 1)/8) + 1):(cnr*(size(subcell_coords, 1)/8))), 1) = sum(temp_reshaped2, 2); 
            end
            %The first size(subcell_coords, 1)/8 control points in
            %subcell_ctrlpts_temp correspond to interpolated control points
            %for all "corner1"s of all the sub-cells, the second set of 
            %size(subcell_coords, 1)/8 control points corresponds to all
            %"corner 2"s, etc. We want to inspect corners 1-8 for each
            %sub-cell, so reshape subcell_ctrlpts_temp so that each column 
            %represents one of the 8 corners and each row represents one 
            %sub-cell.
            subcell_ctrlpts = reshape(subcell_ctrlpts_temp, (size(subcell_coords, 1)/8), 8); 
            %Check interpolated control point signs for each sub-cell (each
            %row of subcell_ctrlpts)
            ctrlpts_signs = sum(sign(subcell_ctrlpts), 2);
            %Check if there are any sub-cells that have all of their
            %control points equal to 0
            zero_ctrlpts_signs = find(~any(sign(subcell_ctrlpts), 2));
            %Extract the row indices of the sub-cells that do NOT have the 
            %same control point signs on all of their corners: these sub-
            %cells will be considered OCCUPIED and their corner coordinates 
            %will be marked for further subdivision in 
            %subcell_coords_occupied, unless lvl_d = b + 1, in which case 
            %the occupied sub-cells will be considered occupied voxels and
            %their corner coordinates will be recorded directly in 
            %reconstructed_vox_pos_corners.
            occupied_subcell_inds = find((abs(ctrlpts_signs) ~= 8));
            occupied_subcell_inds(ismember(occupied_subcell_inds, zero_ctrlpts_signs) == 1) = [];
            if ~isempty(occupied_subcell_inds)
                first_inds = occupied_subcell_inds*8 - 7;
                last_inds = occupied_subcell_inds*8;
                all_inds = [];
                for i = 1:length(first_inds)
                    all_inds((end + 1):(end + 8), 1) = first_inds(i):last_inds(i);
                end
                if lvl_d < b + 1
                    subcell_coords_occupied((end + 1):(end + length(all_inds)), 1:3) = subcell_coords_orig(all_inds, :);
                else
                    reconstructed_vox_pos_corners(rec_vox_pos_cnrs_cntr:(rec_vox_pos_cnrs_cntr + length(all_inds) - 1), 1:3) = subcell_coords_orig(all_inds, :);
                    rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + 8*length(occupied_subcell_inds);
                end
            end %End check if ~isempty(occupied_subcell_inds)
            %If ALL the control points of a sub-cell have the SAME sign, 
            %consider that sub-cell UNoccupied: it will not be subdivided 
            %further, nor recorded as an occupied voxel if lvl_d = b + 1.     
        end %End lvl_d    
    end %End occ_cell
    %profile viewer
end %End lvl

%If there are any surplus rows in reconstructed_vox_pos_corners (because
%the initial matrix size was made too large), remove them now
reconstructed_vox_pos_corners((ismember(reconstructed_vox_pos_corners, [0 0 0], 'rows') == 1), :) = [];

disp('------------------------------------------------------------');
disp('Finding the centre coordinates of all reconstructed occupied voxels ...');
%Our input point cloud's voxel coordinates were considered to be the
%centres of the occupied voxels, so find the midpoint of each voxel in 
%reconstructed_vox_pos_corners: these midpoints will represent our 
%reconstructed voxel positions
reconstructed_vox_pos = zeros((size(reconstructed_vox_pos_corners, 1)/8), 3);
vox_cntr = 1;
for vc = 1:8:(size(reconstructed_vox_pos_corners, 1) - 7)
    %Get the current set of 8 voxel corner coordinates
    vc_coords = reconstructed_vox_pos_corners((vc:(vc + 7)), :);
    %Find the mean of each of the 8 corner coordinates (x, y, and z
    %separately): these mean values represent the centre (x, y, z) location
    %of the current voxel
    reconstructed_vox_pos(vox_cntr, :) = mean(vc_coords, 1);
    vox_cntr = vox_cntr + 1;
end
vox_recon_time = toc;
disp('------------------------------------------------------------');
disp(' ');
disp('************************************************************');
disp(['Time taken to reconstruct voxels: ' num2str(vox_recon_time) ' seconds']);
disp('************************************************************');
disp(' ');

%Plot the reconstructed voxels
figure;
%NOTE: Below, we are only reading in the original PLY file in order 
%to get the corresponding colours assigned to each reconstructed 
%voxel (this will only work when the same number of voxels is 
%reconstructed as in the original point cloud, and these voxels are
%in the same order as the original voxels, which will be done below)
[~, ptcloud, ~] = plyRead(ptcloud_file);
if size(reconstructed_vox_pos, 1) == size(ptcloud, 1)
    %Order the reconstructed voxels according to their Morton codes, so
    %that they are in the same order as the input point cloud at the
    %encoder
    disp('Reordering reconstructed voxels according to Morton codes ...');
    %Get Morton codes for the reconstructed voxel x, y, z coordinates
    mortonCodes = xyzToMorton(reconstructed_vox_pos, b);   %"b" bits for each Morton code
    disp('Morton codes computed');
    %Sort the Morton codes obtained above, in ascending order
    [~, I_vox] = sort(mortonCodes);
    disp('Morton codes sorted');
    %Sort the voxel x, y, z locations in the same order as the sorted 
    %Morton codes
    reconstructed_vox_pos = reconstructed_vox_pos(I_vox, 1:3);
    disp('Reconstructed voxels sorted');
    disp('------------------------------------------------------------');
    %Plot the reconstructed point cloud with the original colours
    %assigned to each reconstructed voxel
    scatter3(reconstructed_vox_pos(:, 1), reconstructed_vox_pos(:, 2), reconstructed_vox_pos(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
else
    %Plot the reconstructed point cloud using a default colour for 
    %all the voxels, since the reconstruction does not contain the same 
    %number of voxels as the original point cloud 
    scatter3(reconstructed_vox_pos(:, 1), reconstructed_vox_pos(:, 2), reconstructed_vox_pos(:, 3), 5, 'filled');
end
axis equal; axis off;
title({'Voxel Reconstruction after using Pruned Octree', 'and Pruned Wavelet Coefficient Tree'});
%Save the above reconstruction as a MATLAB figure and as a PDF image in
%our network directory (NB: The '-bestfit' option maximizes the size of 
%the figure to fill the page, but preserves the aspect ratio of the 
%figure. The figure might not fill the entire page. This option leaves 
%a minimum page margin of .25 inches).
savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_post_pruning']);
print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_post_pruning'], '-dpdf');
disp('Saving reconstructed voxels figure ...');
disp('------------------------------------------------------------');

%For debugging purposes, check if there are any voxels in the original 
%voxelized point cloud that have NOT been reconstructed (i.e., are not 
%present in reconstructed_vox_pos), and if so, plot them
test_vox_diffs = setdiff(ptcloud(:, 1:3), reconstructed_vox_pos, 'rows');
disp(['Total number of original voxels NOT present in the reconstruction: ' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' (' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']);  
if ~isempty(test_vox_diffs)
    figure;
    scatter3(test_vox_diffs(:, 1), test_vox_diffs(:, 2), test_vox_diffs(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
    axis equal; axis off;
    title({'Voxels that were Not Reconstructed,', 'or Not Correctly Reconstructed, at Decoder', ['(' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' = ' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']});
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\missing_voxels_post_pruning']);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\missing_voxels_post_pruning'], '-dpdf');
    disp('Saving missing voxels figure ...');
    disp('------------------------------------------------------------');
end
%Overlay the missing voxels over the reconstructed voxels, to see where
%the gaps are
if ~isempty(test_vox_diffs)
    figure;
    scatter3(reconstructed_vox_pos(:, 1), reconstructed_vox_pos(:, 2), reconstructed_vox_pos(:, 3), 5, 'filled', 'MarkerFaceColor', 'b');
    hold on;
    scatter3(test_vox_diffs(:, 1), test_vox_diffs(:, 2), test_vox_diffs(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
    axis equal; axis off;
    title('Reconstructed and Not-Reconstructed Voxels at Decoder');
    legend('Reconstructed', 'Not Reconstructed', 'Location', 'best');
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\reconstructed_vs_notreconstructed_voxels_post_pruning']);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\reconstructed_vs_notreconstructed_voxels_post_pruning'], '-dpdf');
    disp('Saving reconstructed vs not-reconstructed voxels figure ...');
    disp('------------------------------------------------------------');
end
%For debugging purposes, also check if any voxels are present in 
%reconstructed_vox_pos that were NOT present in the original voxelized 
%point cloud, and if so then plot these
test_vox_diffs2 = setdiff(reconstructed_vox_pos, ptcloud(:, 1:3), 'rows');
disp(['Number of incorrectly reconstructed voxels at decoder: ' num2str(size(test_vox_diffs2, 1))]);
if ~isempty(test_vox_diffs2)
    figure;
    scatter3(test_vox_diffs2(:, 1), test_vox_diffs2(:, 2), test_vox_diffs2(:, 3), 5, 'filled', 'MarkerFaceColor', 'r');
    axis equal; axis off;
    title(['Incorrectly Reconstructed Voxels at Decoder: ' num2str(size(test_vox_diffs2, 1))]);
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_post_pruning']);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_post_pruning'], '-dpdf');
    disp('Saving incorrect voxels figure ...');
    disp('------------------------------------------------------------');
end   
