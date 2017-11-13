function [reconstructed_vox_pos, reconstructed_vox_pos_corners] = voxel_reconstruction_pruning(debug_flag, pp_first_nonempty, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, q_stepsize)

start_vox_recon_time = tic;

%Matrix to store the reconstructed voxel corners: initialize to some large
%value (can shrink later)
reconstructed_vox_pos_corners = zeros(24000000, 3);
%Counter for reconstructed_vox_pos_corners
rec_vox_pos_cnrs_cntr = 1;

%If there are any voxels that have already been reconstructed (because no
%pruning was done further up the octree for these voxels), get their corner
%coordinates
if ~isempty(corner_coords_decoder{b + 1})
    if debug_flag == 1
        disp('Getting corner coordinates for already-reconstructed occupied voxels ...');
        disp(' ');    
    end
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
%decoder-reconstructed control points with the same sign (+/-). Also, store
%the indices and control points of these leaf cells, as we will process 
%them further.
same_sign_leaves = cell(size(pp_leaves_only));
same_sign_leaves_ctrlpts = cell(size(pp_leaves_only));
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
            %Store the index of the current leaf cell
            same_sign_leaves{lvl}((end + 1), 1) = occ_cell;
            %Store the control points of the current leaf cell
            same_sign_leaves_ctrlpts{lvl}(((end + 1):(end + 8)), 1) = current_ctrlpts;
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
    end %End occ_cell
    if debug_flag == 1
        disp(['TOTAL number of leaf cells with all control points having the same sign, at level ' num2str(lvl) ': ' num2str(same_sign_cntr) '/' num2str(length(pp_leaves_only{lvl})) ' (' num2str((same_sign_cntr/(length(pp_leaves_only{lvl}))*100)) '%)']);
        disp(['No. of leaf cells with all 0 control points at level ' num2str(lvl) ': ' num2str(zero_cp_cntr)]);
        disp(' ');
    end
end %End lvl

%Find the first octree level at which same_sign_leaves_ctrlpts is not empty
ssl_first_nonempty = find(~cellfun(@isempty, same_sign_leaves_ctrlpts), 1);
%For all the leaf cells that have all their control points with the same
%sign, check if any of these control points are within +/- q_stepsize/2 of 
%0: if so, change the value of these control points to 0.
for lvl = ssl_first_nonempty:size(same_sign_leaves_ctrlpts, 1)
    if ~isempty(same_sign_leaves_ctrlpts{lvl})
        same_sign_leaves_ctrlpts{lvl}(abs(same_sign_leaves_ctrlpts{lvl}) <= q_stepsize/2) = 0;
    end
end

%For each octree level at which there are leaf cells, except the voxel
%level ...
% figure;
for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    %Temporary: for debugging at one level only
    %parent_indices = [];    %Index of leaf cell at this level, which was used to reconstruct each voxel
    %profile on
    %For each leaf cell at this level ...
    for occ_cell = pp_leaves_only{lvl}'
        if debug_flag == 1
            disp(' ');
            disp(['Reconstructing voxels for leaf cell (occ_cell) ' num2str(occ_cell) ' at level ' num2str(lvl) ':']);
        end
        %Get the 8 corner coordinates of the current leaf cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
%         if (lvl == 10) && (occ_cell == 26121)
%             pause(1);
%         end
%         scatter3(current_corner_coords(:, 1), current_corner_coords(:, 2), current_corner_coords(:, 3), 5, 'filled', 'r');
%         axis equal; axis off;
%         hold on;
        %Get the control points for all 8 corners of the current leaf cell
        if (lvl < ssl_first_nonempty) || ((lvl >= ssl_first_nonempty) && (isempty(find((same_sign_leaves{lvl} == occ_cell), 1))))
            current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        elseif (lvl >= ssl_first_nonempty) && (~isempty(find((same_sign_leaves{lvl} == occ_cell), 1)))
            %Get the index of the current occ_cell inside same_sign_leaves
            ind = find((same_sign_leaves{lvl} == occ_cell), 1);
            current_ctrlpts = same_sign_leaves_ctrlpts{lvl}((ind*8 - 7):(ind*8));
        end
        %If current_ctrlpts are all 0, consider every voxel that belongs to
        %the current leaf cell as being occupied, and do not process the
        %current leaf cell any further
        if ~any(current_ctrlpts)
            %Subdivide the current occ_cell into 1 x 1 x 1 voxels
            min_x = min(current_corner_coords(:, 1));  
            min_y = min(current_corner_coords(:, 2));
            min_z = min(current_corner_coords(:, 3));
            max_x = max(current_corner_coords(:, 1));
            max_y = max(current_corner_coords(:, 2));
            max_z = max(current_corner_coords(:, 3));
            %Offsets for all 8 corners of each voxel
            offsets = [0 0 0;
                1 0 0;
                1 1 0;
                0 1 0;
                0 0 1;
                1 0 1;
                1 1 1;
                0 1 1]; 
            %For each shift by 1 on the x-axis
            for x_ax = min_x:(max_x - 1)
                %For each shift by 1 on the y-axis
                for y_ax = min_y:(max_y - 1)
                    %For each shift by 1 on the z-axis
                    for z_ax = min_z:(max_z - 1)
                        %Compute all 8 corners of the voxel whose first
                        %corner is at (x_ax, y_ax, z_ax)
                        reconstructed_vox_pos_corners(rec_vox_pos_cnrs_cntr:(rec_vox_pos_cnrs_cntr + 8 - 1), 1:3) = [x_ax y_ax z_ax] + offsets;      
                        rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + 8;
                    end
                end %End y_ax
            end %End x_ax
            %Do not process the current occ_cell any further
            continue;
        end 
        %Keep subdividing the current leaf cell until we reach cells of 
        %size 1 x 1 x 1 (i.e., voxels). For each sub-cell at each level,
        %interpolate between the leaf cell's (occ_cell's) control points.
        %The sub-cells that have interpolated control points with different 
        %signs will be subdivided further; at the end, voxels that have 
        %interpolated control points with different signs (or voxels that
        %have all 0 control points) will be considered occupied.
        for lvl_d = (lvl + 1):(b + 1)
            %Compute the corner coordinates of each of the sub-cells at
            %lvl_d ...
            if lvl_d == (lvl + 1)
                if debug_flag == 1
                    disp(['--Processing sub-cells at lvl_d ' num2str(lvl_d) ': 8 sub-cells to process at current lvl_d']);
                end
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
                if debug_flag == 1
                    disp(['--Processing sub-cells at lvl_d ' num2str(lvl_d) ': ' num2str(size(subcell_coords_occupied, 1)/8) ' occupied sub-cells at previous lvl_d, so ' num2str(size(subcell_coords_occupied, 1)) ' sub-cells to process at current lvl_d']);    
                end
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
            if debug_flag == 1
                disp(['  Computing corner coordinates for each sub-cell at this level (' num2str(size(zc_coords, 1)) ' sub-cells in total)']);           
            end
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
            if debug_flag == 1
                disp(['  Computing interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
            end
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
%             %If any subcell_ctrlpts are small (within +/- q_stepsize/2),
%             %set them to be 0
%             subcell_ctrlpts(abs(subcell_ctrlpts) <= q_stepsize/2) = 0;
            %Check if any of the current sub-cells have all 0 control
            %points: in this case, consider every voxel that belongs to
            %that sub-cell as being occupied, and do not process this
            %sub-cell any further
            all_zero_subcells = find(~any(subcell_ctrlpts, 2));
            if ~isempty(all_zero_subcells)
                for az_sc = all_zero_subcells'
                    %For the current all-zero sub-cell, find its corner
                    %coordinates
                    cnr_coords = subcell_coords_orig(((az_sc*8 - 7):(az_sc*8)), :);
                    if lvl_d < b + 1
                        %Subdivide the current sub-cell into 1 x 1 x 1 
                        %voxels
                        sc_min_x = min(cnr_coords(:, 1));  
                        sc_min_y = min(cnr_coords(:, 2));
                        sc_min_z = min(cnr_coords(:, 3));
                        sc_max_x = max(cnr_coords(:, 1));
                        sc_max_y = max(cnr_coords(:, 2));
                        sc_max_z = max(cnr_coords(:, 3));
                        %Offsets for all 8 corners of each voxel
                        offsets = [0 0 0;
                            1 0 0;
                            1 1 0;
                            0 1 0;
                            0 0 1;
                            1 0 1;
                            1 1 1;
                            0 1 1]; 
                        %For each shift by 1 on the x-axis
                        for x_ax = sc_min_x:(sc_max_x - 1)
                            %For each shift by 1 on the y-axis
                            for y_ax = sc_min_y:(sc_max_y - 1)
                                %For each shift by 1 on the z-axis
                                for z_ax = sc_min_z:(sc_max_z - 1)
                                    %Compute all 8 corners of the voxel 
                                    %whose first corner is at 
                                    %(x_ax, y_ax, z_ax)
                                    reconstructed_vox_pos_corners(rec_vox_pos_cnrs_cntr:(rec_vox_pos_cnrs_cntr + 8 - 1), 1:3) = [x_ax y_ax z_ax] + offsets;      
                                    rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + 8;
                                end
                            end %End y_ax
                        end %End x_ax
                    else
                        %Record the corners of the current voxel, as this
                        %voxel is considered to be occupied
                        reconstructed_vox_pos_corners(rec_vox_pos_cnrs_cntr:(rec_vox_pos_cnrs_cntr + 8 - 1), 1:3) = cnr_coords;      
                        rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + 8;
                    end   
                end %End az_sc
            end %End check if ~isempty(all_zero_subcells)
            %Check interpolated control point signs for each sub-cell (each
            %row of subcell_ctrlpts)
            ctrlpts_signs = sum(sign(subcell_ctrlpts), 2);
            %Extract the row indices of the sub-cells that do NOT have the 
            %same control point signs on all of their corners: these sub-
            %cells will be considered OCCUPIED and their corner coordinates 
            %will be marked for further subdivision in 
            %subcell_coords_occupied, unless lvl_d = b + 1, in which case 
            %the occupied sub-cells will be considered occupied voxels and
            %their corner coordinates will be recorded directly in 
            %reconstructed_vox_pos_corners.
            occupied_subcell_inds = find((abs(ctrlpts_signs) ~= 8));
            %Exclude the indices of any sub-cells that have all zero
            %control points, as these have already been processed above
            occupied_subcell_inds(ismember(occupied_subcell_inds, all_zero_subcells) == 1) = [];
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
                    %parent_indices(rec_vox_pos_cnrs_cntr:(rec_vox_pos_cnrs_cntr + length(all_inds) - 1), 1) = repmat(occ_cell, size(subcell_coords_orig(all_inds, :), 1), 1); %For debugging only
                    rec_vox_pos_cnrs_cntr = rec_vox_pos_cnrs_cntr + length(all_inds);
                end
            end %End check if ~isempty(occupied_subcell_inds)
            %If ALL the control points of a sub-cell have the SAME sign, 
            %consider that sub-cell UNoccupied: it will not be subdivided 
            %further, nor recorded as an occupied voxel if lvl_d = b + 1.     
        end %End lvl_d   
%         if occ_cell >= 102
%             %test = reconstructed_vox_pos_corners((1:(rec_vox_pos_cnrs_cntr - 1)), :);
%             test = subcell_coords_orig(all_inds, :);
%             test_means = zeros((size(test, 1)/8), 3);
%             vox_cntr = 1;
%             for vc = 1:8:(size(test, 1) - 7)
%                 %Get the current set of 8 voxel corner coordinates
%                 vc_coords = test((vc:(vc + 7)), :);
%                 %Find the mean of each of the 8 corner coordinates (x, y, and z
%                 %separately): these mean values represent the centre (x, y, z) location
%                 %of the current voxel
%                 test_means(vox_cntr, :) = mean(vc_coords, 1);
%                 vox_cntr = vox_cntr + 1;
%             end
%             scatter3(test_means(:, 1), test_means(:, 2), test_means(:, 3), 5, 'filled', 'b');
%             hold on;
%             scatter3(current_corner_coords(:, 1), current_corner_coords(:, 2), current_corner_coords(:, 3), 5, 'filled', 'm');
%             axis equal; axis off;
%         end
    end %End occ_cell
    %profile viewer
end %End lvl

%If there are any surplus rows in reconstructed_vox_pos_corners (because
%the initial matrix size was made too large), remove them now
%reconstructed_vox_pos_corners((ismember(reconstructed_vox_pos_corners, [0 0 0], 'rows') == 1), :) = [];
reconstructed_vox_pos_corners((rec_vox_pos_cnrs_cntr:end), :) = [];

if debug_flag == 1
    disp('------------------------------------------------------------');
    disp('Finding the centre coordinates (midpoints) of all reconstructed occupied voxels ...');
end
%Our input point cloud's voxel coordinates were considered to be the
%centres of the occupied voxels, so find the midpoint of each voxel in 
%reconstructed_vox_pos_corners: these midpoints will represent our 
%reconstructed voxel positions
reconstructed_vox_pos = zeros((size(reconstructed_vox_pos_corners, 1)/8), 3);
%parent_indices_final = [];  %For debugging only: leaf cell indices for reconstructed voxels
vox_cntr = 1;
for vc = 1:8:(size(reconstructed_vox_pos_corners, 1) - 7)
    %Get the current set of 8 voxel corner coordinates
    vc_coords = reconstructed_vox_pos_corners((vc:(vc + 7)), :);
    %Find the mean of each of the 8 corner coordinates (x, y, and z
    %separately): these mean values represent the centre (x, y, z) location
    %of the current voxel
    reconstructed_vox_pos(vox_cntr, :) = mean(vc_coords, 1);
    %parent_indices_final(vox_cntr, 1) = parent_indices(vc);
    vox_cntr = vox_cntr + 1;
end
vox_recon_time = toc(start_vox_recon_time);
if debug_flag == 1
    disp('------------------------------------------------------------');
end
disp(' ');
disp('************************************************************');
disp(['Time taken to reconstruct voxels: ' num2str(vox_recon_time) ' seconds']);
disp('************************************************************');
