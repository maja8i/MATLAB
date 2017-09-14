function [reconstructed_vox_pos, reconstructed_vox_pos_corners] = voxel_reconstruction_pruning(pp_first_nonempty, SpatialIndex, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, ptcloud_file, ptcloud_name)

%Matrix to store the reconstructed voxel corners
reconstructed_vox_pos_corners = [];

%If there are any voxels that have already been reconstructed (because no
%pruning was done further up the octree for these voxels), get their corner
%coordinates
if ~isempty(corner_coords_decoder{b + 1})
    disp(' ');
    disp('Getting corner coordinates for already-reconstructed occupied voxels ...');
    disp(' ');    
    reconstructed_vox_pos_corners = corner_coords_decoder{b + 1};
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

%For each octree level at which there are leaf cells
for lvl = pp_first_nonempty:size(SpatialIndex, 1)
    %For each leaf cell at this level
    for occ_cell = pp_leaves_only{lvl}'
        disp(['Reconstructing voxels for leaf cell (occ_cell) ' num2str(occ_cell) ' at level ' num2str(lvl) ':']);
        %Get the 8 corner coordinates of the current leaf cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
        %Get the control points for all 8 corners of the current leaf cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Keep subdividing the current leaf cell until we reach cells of 
        %size 1x1x1 (i.e., voxels). For each sub-cell at each level,
        %interpolate between the leaf cell's (occ_cell's) control points.
        %The voxels that have interpolated control points with different 
        %signs at the end will be considered occupied.
        for lvl_d = (lvl + 1):(b + 1)
            %profile on
            %Compute the corner coordinates of each of the 8 sub-cells at
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
            else
                disp(['--Processing sub-cells at lvl_d ' num2str(lvl_d) ': ' num2str(size(subcell_coords_occupied, 1)/8) ' occupied sub-cells at previous lvl_d, so ' num2str(size(subcell_coords_occupied, 1)) ' sub-cells to process at current lvl_d']);    
                %Find the minimum, maximum and midpoint x, y, z of only the
                %cells at the previous level that have been marked as
                %candidates for further subdivision (those stored in 
                %subcell_coords_occupied)
                min_x = zeros((size(subcell_coords_occupied, 1))/8, 1);
                min_y = zeros((size(subcell_coords_occupied, 1))/8, 1);
                min_z = zeros((size(subcell_coords_occupied, 1))/8, 1);
                max_x = zeros((size(subcell_coords_occupied, 1))/8, 1);
                max_y = zeros((size(subcell_coords_occupied, 1))/8, 1);
                max_z = zeros((size(subcell_coords_occupied, 1))/8, 1);
                mid_x = zeros((size(subcell_coords_occupied, 1))/8, 1);
                mid_y = zeros((size(subcell_coords_occupied, 1))/8, 1);
                mid_z = zeros((size(subcell_coords_occupied, 1))/8, 1);
                %There will be one min, one max and one mid point (for x, 
                %y, z separately) per 8 sub-cells at the current lvl_d 
                %(there are 8 sub-cells, not all necessarily occupied, per 
                %occupied cell at subcell_coords_occupied)
                cell_cntr = 1;
                for sub_cell_8set = 1:8:((size(subcell_coords_occupied, 1)) - 7)
                    min_x(cell_cntr) = min(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 1));
                    min_y(cell_cntr) = min(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 2));
                    min_z(cell_cntr) = min(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 3));
                    max_x(cell_cntr) = max(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 1));
                    max_y(cell_cntr) = max(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 2));
                    max_z(cell_cntr) = max(subcell_coords_occupied((sub_cell_8set:(sub_cell_8set + 7)), 3));
                    mid_x(cell_cntr) = (min_x(cell_cntr) + max_x(cell_cntr))/2;
                    mid_y(cell_cntr) = (min_y(cell_cntr) + max_y(cell_cntr))/2;
                    mid_z(cell_cntr) = (min_z(cell_cntr) + max_z(cell_cntr))/2;
                    cell_cntr = cell_cntr + 1;
                end  
            end %End check if lvl_d == (lvl + 1) 
            %Use the min, max, and mid points found above, to compute the
            %corner coordinates of each of the 8 sub-cells resulting from
            %subdividing each of the cells whose coordinates are stored in
            %zc_coords, below
            if (lvl_d == lvl + 1)
                zc_coords = current_corner_coords;  %"zc" stands for zero crossing, as we only subdivide occupied cells
            elseif (lvl_d > lvl + 1)
                zc_coords = subcell_coords_occupied;
            end
            subcell_coords = zeros(size(zc_coords, 1)*8, 3); %Temporary matrix to store the corner coordinates of each sub-cell at the current level 
            scc_cntr = 1;   %Counter for corner coordinates of the sub-cells at the current level
            cell_cntr = 1;
            disp(['  Computing corner coordinates for each sub-cell at this level (' num2str(size(zc_coords, 1)) ' sub-cells in total)']);
            for sub_cell_8set = 1:8:(size(zc_coords, 1) - 7)
                %Sub-cell 1
                if lvl_d == (lvl + 1)
                    subcell_coords((scc_cntr:(scc_cntr + 7)), :) = [min_x min_y min_z; mid_x min_y min_z; mid_x mid_y min_z; min_x mid_y min_z; min_x min_y mid_z; mid_x min_y mid_z; mid_x mid_y mid_z; min_x mid_y mid_z];   
                else
                    subcell_coords((scc_cntr:(scc_cntr + 7)), :) = [min_x(cell_cntr) min_y(cell_cntr) min_z(cell_cntr); mid_x(cell_cntr) min_y(cell_cntr) min_z(cell_cntr); mid_x(cell_cntr) mid_y(cell_cntr) min_z(cell_cntr); min_x(cell_cntr) mid_y(cell_cntr) min_z(cell_cntr); min_x(cell_cntr) min_y(cell_cntr) mid_z(cell_cntr); mid_x(cell_cntr) min_y(cell_cntr) mid_z(cell_cntr); mid_x(cell_cntr) mid_y(cell_cntr) mid_z(cell_cntr); min_x(cell_cntr) mid_y(cell_cntr) mid_z(cell_cntr)];   
                end
                %Sub-cell 2
                subcell_coords(((scc_cntr + 8):(scc_cntr + 15)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [2^(b + 1 - lvl_d) 0 0];
                %Sub-cell 3
                subcell_coords(((scc_cntr + 16):(scc_cntr + 23)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d) 0];
                %Sub-cell 4
                subcell_coords(((scc_cntr + 24):(scc_cntr + 31)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [0 2^(b + 1 - lvl_d) 0];
                %Sub-cell 5
                subcell_coords(((scc_cntr + 32):(scc_cntr + 39)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [0 0 2^(b + 1 - lvl_d)];
                %Sub-cell 6
                subcell_coords(((scc_cntr + 40):(scc_cntr + 47)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [2^(b + 1 - lvl_d) 0 2^(b + 1 - lvl_d)];
                %Sub-cell 7
                subcell_coords(((scc_cntr + 48):(scc_cntr + 55)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d)];
                %Sub-cell 8
                subcell_coords(((scc_cntr + 56):(scc_cntr + 63)), :) = subcell_coords((scc_cntr:(scc_cntr + 7)), :) + [0 2^(b + 1 - lvl_d) 2^(b + 1 - lvl_d)];
                %Increment scc_cntr for the next set of 8 sub-cells
                scc_cntr = scc_cntr + 64; 
                cell_cntr = cell_cntr + 1;
            end %End sub_cell_8set
            %Array that will store the corner coordinates of subcells at
            %the current lvl_d, which are candidates for further
            %subdivision (if lvl_d < b + 1)
            subcell_coords_occupied = [];
            %For each corner computed in subcell_coords, calculate the
            %interpolated control point for this corner, by trilinearly
            %interpolating between the original parent control points,
            %current_ctrlpts
            disp(['  Computing interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
            for sub_cell_8set = 1:8:(size(subcell_coords, 1) - 7)
                %Get all 8 corner coordinates for the current sub-cell
                sc_coords = subcell_coords((sub_cell_8set:(sub_cell_8set + 7)), :);
                %Normalize sc_coords to be in the range [0, 1], because the
                %trilinear interpolation formula (used below) expects the 
                %x, y, z values to be in this range
                sc_coords_orig = sc_coords;
                sc_coords = (sc_coords - current_corner_coords(1, :))/(2^(b + 1 - lvl));    %Subtract the origin of the original parent cell, and divide by the cell width
                %Initialize an array to store the interpolated control
                %points for all 8 corners of the current sub-cell
                subcell_ctrlpts = zeros(8, 1);
                %Compute all 8 control points for the corners of the 
                %current sub-cell
                mult_matrix = [(1 - sc_coords(1, 1))*(1 - sc_coords(1, 2))*(1 - sc_coords(1, 3));
                    sc_coords(1, 1)*(1 - sc_coords(1, 2))*(1 - sc_coords(1, 3));
                    sc_coords(1, 1)*sc_coords(1, 2)*(1 - sc_coords(1, 3));
                    (1 - sc_coords(1, 1))*sc_coords(1, 2)*(1 - sc_coords(1, 3));
                    (1 - sc_coords(1, 1))*(1 - sc_coords(1, 2))*sc_coords(1, 3);
                    sc_coords(1, 1)*(1 - sc_coords(1, 2))*sc_coords(1, 3);
                    sc_coords(1, 1)*sc_coords(1, 2)*sc_coords(1, 3);
                    (1 - sc_coords(1, 1))*sc_coords(1, 2)*sc_coords(1, 3);  %End corner 1
                    (1 - sc_coords(2, 1))*(1 - sc_coords(2, 2))*(1 - sc_coords(2, 3));
                    sc_coords(2, 1)*(1 - sc_coords(2, 2))*(1 - sc_coords(2, 3));
                    sc_coords(2, 1)*sc_coords(2, 2)*(1 - sc_coords(2, 3));
                    (1 - sc_coords(2, 1))*sc_coords(2, 2)*(1 - sc_coords(2, 3));
                    (1 - sc_coords(2, 1))*(1 - sc_coords(2, 2))*sc_coords(2, 3);
                    sc_coords(2, 1)*(1 - sc_coords(2, 2))*sc_coords(2, 3);
                    sc_coords(2, 1)*sc_coords(2, 2)*sc_coords(2, 3);
                    (1 - sc_coords(2, 1))*sc_coords(2, 2)*sc_coords(2, 3);  %End corner 2
                    (1 - sc_coords(3, 1))*(1 - sc_coords(3, 2))*(1 - sc_coords(3, 3));
                    sc_coords(3, 1)*(1 - sc_coords(3, 2))*(1 - sc_coords(3, 3));
                    sc_coords(3, 1)*sc_coords(3, 2)*(1 - sc_coords(3, 3));
                    (1 - sc_coords(3, 1))*sc_coords(3, 2)*(1 - sc_coords(3, 3));
                    (1 - sc_coords(3, 1))*(1 - sc_coords(3, 2))*sc_coords(3, 3);
                    sc_coords(3, 1)*(1 - sc_coords(3, 2))*sc_coords(3, 3);
                    sc_coords(3, 1)*sc_coords(3, 2)*sc_coords(3, 3);
                    (1 - sc_coords(3, 1))*sc_coords(3, 2)*sc_coords(3, 3);  %End corner 3
                    (1 - sc_coords(4, 1))*(1 - sc_coords(4, 2))*(1 - sc_coords(4, 3));
                    sc_coords(4, 1)*(1 - sc_coords(4, 2))*(1 - sc_coords(4, 3));
                    sc_coords(4, 1)*sc_coords(4, 2)*(1 - sc_coords(4, 3));
                    (1 - sc_coords(4, 1))*sc_coords(4, 2)*(1 - sc_coords(4, 3));
                    (1 - sc_coords(4, 1))*(1 - sc_coords(4, 2))*sc_coords(4, 3);
                    sc_coords(4, 1)*(1 - sc_coords(4, 2))*sc_coords(4, 3);
                    sc_coords(4, 1)*sc_coords(4, 2)*sc_coords(4, 3);
                    (1 - sc_coords(4, 1))*sc_coords(4, 2)*sc_coords(4, 3);  %End corner 4
                    (1 - sc_coords(5, 1))*(1 - sc_coords(5, 2))*(1 - sc_coords(5, 3));
                    sc_coords(5, 1)*(1 - sc_coords(5, 2))*(1 - sc_coords(5, 3));
                    sc_coords(5, 1)*sc_coords(5, 2)*(1 - sc_coords(5, 3));
                    (1 - sc_coords(5, 1))*sc_coords(5, 2)*(1 - sc_coords(5, 3));
                    (1 - sc_coords(5, 1))*(1 - sc_coords(5, 2))*sc_coords(5, 3);
                    sc_coords(5, 1)*(1 - sc_coords(5, 2))*sc_coords(5, 3);
                    sc_coords(5, 1)*sc_coords(5, 2)*sc_coords(5, 3);
                    (1 - sc_coords(5, 1))*sc_coords(5, 2)*sc_coords(5, 3);  %End corner 5
                    (1 - sc_coords(6, 1))*(1 - sc_coords(6, 2))*(1 - sc_coords(6, 3));
                    sc_coords(6, 1)*(1 - sc_coords(6, 2))*(1 - sc_coords(6, 3));
                    sc_coords(6, 1)*sc_coords(6, 2)*(1 - sc_coords(6, 3));
                    (1 - sc_coords(6, 1))*sc_coords(6, 2)*(1 - sc_coords(6, 3));
                    (1 - sc_coords(6, 1))*(1 - sc_coords(6, 2))*sc_coords(6, 3);
                    sc_coords(6, 1)*(1 - sc_coords(6, 2))*sc_coords(6, 3);
                    sc_coords(6, 1)*sc_coords(6, 2)*sc_coords(6, 3);
                    (1 - sc_coords(6, 1))*sc_coords(6, 2)*sc_coords(6, 3);  %End corner 6
                    (1 - sc_coords(7, 1))*(1 - sc_coords(7, 2))*(1 - sc_coords(7, 3));
                    sc_coords(7, 1)*(1 - sc_coords(7, 2))*(1 - sc_coords(7, 3));
                    sc_coords(7, 1)*sc_coords(7, 2)*(1 - sc_coords(7, 3));
                    (1 - sc_coords(7, 1))*sc_coords(7, 2)*(1 - sc_coords(7, 3));
                    (1 - sc_coords(7, 1))*(1 - sc_coords(7, 2))*sc_coords(7, 3);
                    sc_coords(7, 1)*(1 - sc_coords(7, 2))*sc_coords(7, 3);
                    sc_coords(7, 1)*sc_coords(7, 2)*sc_coords(7, 3);
                    (1 - sc_coords(7, 1))*sc_coords(7, 2)*sc_coords(7, 3);  %End corner 7
                    (1 - sc_coords(8, 1))*(1 - sc_coords(8, 2))*(1 - sc_coords(8, 3));
                    sc_coords(8, 1)*(1 - sc_coords(8, 2))*(1 - sc_coords(8, 3));
                    sc_coords(8, 1)*sc_coords(8, 2)*(1 - sc_coords(8, 3));
                    (1 - sc_coords(8, 1))*sc_coords(8, 2)*(1 - sc_coords(8, 3));
                    (1 - sc_coords(8, 1))*(1 - sc_coords(8, 2))*sc_coords(8, 3);
                    sc_coords(8, 1)*(1 - sc_coords(8, 2))*sc_coords(8, 3);
                    sc_coords(8, 1)*sc_coords(8, 2)*sc_coords(8, 3);
                    (1 - sc_coords(8, 1))*sc_coords(8, 2)*sc_coords(8, 3)]; %End corner 8
                temp = repmat(current_ctrlpts, 8, 1).*mult_matrix;
                subcell_ctrlpts_cntr = 1;
                for c = 1:8:(size(temp, 1) - 7)
                    subcell_ctrlpts(subcell_ctrlpts_cntr) = sum(temp(c:(c + 7)));
                    subcell_ctrlpts_cntr = subcell_ctrlpts_cntr + 1;
                end
                subcell_ctrlpts = subcell_ctrlpts'; %Want a column vector
                %Check interpolated control point signs
                if (sum(sign(subcell_ctrlpts)) ~= length(subcell_ctrlpts))&&(sum(sign(subcell_ctrlpts)) ~= -length(subcell_ctrlpts))
                    %If the control points do NOT all have the same sign, 
                    %consider the current sub-cell OCCUPIED: it will be 
                    %marked as a candidate for further subdivision in
                    %subcell_coords_occupied, unless lvl_d = b + 1, in 
                    %which case it will be considered an occupied voxel and 
                    %recorded in reconstructed_vox_pos_corners
                    if lvl_d < b + 1
                        subcell_coords_occupied(((end + 1):(end + 8)), 1:3) = sc_coords_orig;
                    else
                        reconstructed_vox_pos_corners(((end + 1):(end + 8)), 1:3) = sc_coords_orig;
                    end
                end
                %If all the control points have the SAME sign, consider the
                %current sub-cell UNoccupied: it will not be subdivided 
                %further, nor recorded as an occupied voxel if 
                %lvl_d = b + 1. Continue with checking the next sub-cell.
            end %End sub_cell_8set        
            %profile viewer
        end %End lvl_d    
    end %End occ_cell
end %End lvl

disp('------------------------------------------------------------');
disp('Finding the centre coordinates of all reconstructed occupied voxels ...');
%Our input point cloud's voxel coordinates were considered to be the
%centres of the 1x1x1 voxels, so find the midpoint of each voxel in 
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
disp('------------------------------------------------------------');

%Plot the reconstructed voxels
figure;
%NOTE: Below, we are only reading in the original PLY file in order 
%to get the corresponding colours assigned to each reconstructed 
%voxel (this will only work when the same number of voxels is 
%reconstructed as in the original point cloud, and these voxels are
%in the same order as the original voxels, which was done above)
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
    %all voxels, since the reconstruction does not contain the same 
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
%voxelized point cloud that have not been reconstructed (i.e., are not 
%present in reconstructed_vox_pos), and if so, plot them
test_vox_diffs = setdiff(ptcloud(:, 1:3), reconstructed_vox_pos, 'rows');
disp(['Total number of missing voxels or incorrectly reconstructed voxels in reconstruction at decoder: ' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' (' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']);  
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
