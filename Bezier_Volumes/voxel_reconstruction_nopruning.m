function [reconstructed_vox_pos, reconstructed_vox_pos_corners, subcell_coords_all] = voxel_reconstruction_nopruning(SpatialIndex, corner_coords_decoder, reconstruction_decoder, ctrl_pts_pointers, b, max_lvl, ptcloud_file, ptcloud_name)

%For a chosen (set of) octree level(s), reconstruct the voxel positions 
%that define the surface (shape) of the input 3D point cloud, by using
%only the reconstructed control points at the chosen octree level(s) 
%and interpolating between them for the higher octree levels, to figure 
%out through which octree cells the surface passes, down to the leaf 
%level. Finally, at the leaf level, the midpoints of the voxels that 
%contain the zero crossings will be the reconstructed voxel positions.
%for lvl = [4, 5, 6, 7, 8]
for lvl = max_lvl
    disp(['Reconstructing voxels for level ' num2str(lvl) ' control points ...']);
    disp(' ');
    %Cell array to store the coordinates of octree cells at all levels  
    %from (lvl + 1) to the leaves, which are candidates for further 
    %subdivision or, at the leaf level, contain zero crossings
    subcell_coords_all = cell((b + 1), size(SpatialIndex{lvl}, 1));  
    %For each occupied cell at the current octree level, lvl ...
    for occ_cell = 1:size(SpatialIndex{lvl}, 1) %We can only use this for levels for which we know how many occupied cells there are
        disp(['Reconstructing voxels for occ_cell ' num2str(occ_cell) '/' num2str(size(SpatialIndex{lvl}, 1)) ':']);
        %Get the control points for all 8 corners of this cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Get the corner coordinates for each of the 8 corners of this 
        %cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
        %If we are at the leaf level ...
        if lvl == (b + 1)
            %Just get the current occupied voxel's corner coordinates 
            subcell_coords_all{lvl, occ_cell} = current_corner_coords;
            continue;
        %If we are NOT at the leaf level ...
        else 
            %For each descendant of the current occ_cell, right down to
            %the leaves ...
            for lvl_d = (lvl + 1):(b + 1)     
                disp(['--Processing sub-cells at level ' num2str(lvl_d)]);
                %Find the unique minimum, unique maximum, and midpoint
                %coordinates for the x, y, and z dimensions of the 
                %sub-cells at this level
                if lvl_d == (lvl + 1)
                    unq_min_x = min(unique(current_corner_coords(:, 1)));
                    unq_min_y = min(unique(current_corner_coords(:, 2)));
                    unq_min_z = min(unique(current_corner_coords(:, 3)));
                    unq_max_x = max(unique(current_corner_coords(:, 1)));
                    unq_max_y = max(unique(current_corner_coords(:, 2)));
                    unq_max_z = max(unique(current_corner_coords(:, 3)));
                    mid_x = (unq_min_x + unq_max_x)/2;
                    mid_y = (unq_min_y + unq_max_y)/2;
                    mid_z = (unq_min_z + unq_max_z)/2;
                else
                    %There will be one unq_min, one unq_max, and one 
                    %mid point (for x, y, z separately) per 8 sub-cells 
                    %at this level (there are 8 sub-cells (not all
                    %necessarily occupied) per occupied cell at the 
                    %previous level). The coordinates of cells at the 
                    %previous level, which are candidates for further 
                    %subdivision, are stored in 
                    %subcell_coords_all(lvl_d - 1).
                    unq_min_x = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    unq_min_y = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    unq_min_z = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    unq_max_x = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    unq_max_y = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    unq_max_z = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    mid_x = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    mid_y = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    mid_z = zeros((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1))/8, 1);
                    cell_cntr = 1;
                    for sub_cell_8set = 1:8:((size(subcell_coords_all{(lvl_d - 1), occ_cell}, 1)) - 7)
                        %unq_min_x(cell_cntr) = min(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 1))) + 2^(b + 1 - lvl_d)*(sub_cell_8set - 1);
                        unq_min_x(cell_cntr) = min(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 1)));
                        unq_min_y(cell_cntr) = min(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 2)));
                        unq_min_z(cell_cntr) = min(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 3)));
                        unq_max_x(cell_cntr) = max(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 1)));
                        unq_max_y(cell_cntr) = max(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 2)));
                        unq_max_z(cell_cntr) = max(unique(subcell_coords_all{(lvl_d - 1), occ_cell}((sub_cell_8set:(sub_cell_8set + 7)), 3)));
                        mid_x(cell_cntr) = (unq_min_x(cell_cntr) + unq_max_x(cell_cntr))/2;
                        mid_y(cell_cntr) = (unq_min_y(cell_cntr) + unq_max_y(cell_cntr))/2;
                        mid_z(cell_cntr) = (unq_min_z(cell_cntr) + unq_max_z(cell_cntr))/2;
                        cell_cntr = cell_cntr + 1;
                    end
                end %End check if lvl_d == (lvl + 1)                  
                %Compute the 8 corner coordinates for each of the 8 
                %sub-cells resulting from subdividing each of the cells 
                %whose coordinates are stored in subcell_coords_all 
                %(if lvl_d > lvl + 1) or in current_corner_coords (if 
                %lvl_d = lvl + 1)
                if (lvl_d == lvl + 1)
                    zc_coords = current_corner_coords;  %"zc" stands for zero crossing
                elseif (lvl_d > lvl + 1)
                    zc_coords = subcell_coords_all{(lvl_d - 1), occ_cell};
                end
                subcell_coords = zeros(size(zc_coords, 1)*8, 3); %Matrix to store the corner coordinates of each sub-cell at the current level 
                scc_cntr = 1;   %Counter for corner coordinates of the sub-cells at the current level
                cell_cntr = 1;
                disp(['  Computing corner coordinates for each sub-cell at this level (' num2str(size(zc_coords, 1)) ' sub-cells in total)']);
                for sub_cell_8set = 1:8:(size(zc_coords, 1) - 7)
                    %Sub-cell 1
                    if lvl_d == (lvl + 1)
                        subcell_coords((scc_cntr:(scc_cntr + 7)), :) = [unq_min_x unq_min_y unq_min_z; mid_x unq_min_y unq_min_z; mid_x mid_y unq_min_z; unq_min_x mid_y unq_min_z; unq_min_x unq_min_y mid_z; mid_x unq_min_y mid_z; mid_x mid_y mid_z; unq_min_x mid_y mid_z];   
                    else
                        subcell_coords((scc_cntr:(scc_cntr + 7)), :) = [unq_min_x(cell_cntr) unq_min_y(cell_cntr) unq_min_z(cell_cntr); mid_x(cell_cntr) unq_min_y(cell_cntr) unq_min_z(cell_cntr); mid_x(cell_cntr) mid_y(cell_cntr) unq_min_z(cell_cntr); unq_min_x(cell_cntr) mid_y(cell_cntr) unq_min_z(cell_cntr); unq_min_x(cell_cntr) unq_min_y(cell_cntr) mid_z(cell_cntr); mid_x(cell_cntr) unq_min_y(cell_cntr) mid_z(cell_cntr); mid_x(cell_cntr) mid_y(cell_cntr) mid_z(cell_cntr); unq_min_x(cell_cntr) mid_y(cell_cntr) mid_z(cell_cntr)];   
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
                %For each corner coordinate in subcell_coords, compute
                %the control point associated with this corner, by 
                %interpolating (tri-linear interpolation) between the 
                %original parent control points, current_ctrlpts
                disp(['  Computing interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
                zcc_coords_cntr = 1;    %Counter for coordinates of cells at lvl_d, which will be candidates for further subdivision
                for sub_cell_8set = 1:8:(size(subcell_coords, 1) - 7)
                    %Get all 8 corner coordinates for the current sub-cell
                    sc_coords = subcell_coords((sub_cell_8set:(sub_cell_8set + 7)), :);
                    %Normalize sc_coords to be in the range [0, 1],  
                    %because the trilinear interpolation formula (used 
                    %below) expects the (x, y, z) values to be in this 
                    %range
                    sc_coords_orig = sc_coords;
                    %sc_coords = (sc_coords - sc_coords(1, :))/(2^(b + 1 - lvl_d));   %Subtract the origin of the current sub-cell, and divide by the cell width
                    sc_coords = (sc_coords - current_corner_coords(1, :))/(2^(b + 1 - lvl));
                    %Initialize an array to store the interpolated 
                    %control points for all 8 corners of the current 
                    %sub-cell
                    subcell_ctrlpts = zeros(8, 1);
                    %Compute all 8 control points for the corners of  
                    %the current sub-cell
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

                    %If we are one level before the leaf level, or at
                    %the leaf level
                    if lvl_d >= b
                        %Check interpolated control point signs
                        if (sum(sign(subcell_ctrlpts)) == length(subcell_ctrlpts))||(sum(sign(subcell_ctrlpts)) == -length(subcell_ctrlpts))
                            %If all the control points have the same 
                            %sign, consider the current sub-cell 
                            %unoccupied: it will not be subdivided 
                            %further
                            continue;
                        else
                            %If the control points do not all have the 
                            %same sign, consider the current sub-cell 
                            %occupied: it will be subdivided further
                            %unless we are at the leaf level (b + 1)
                            subcell_coords_all{lvl_d, occ_cell}((zcc_coords_cntr:(zcc_coords_cntr + 7)), 1:3) = sc_coords_orig;
                            zcc_coords_cntr = zcc_coords_cntr + 8;
                        end
                    %If we are not at the leaf level or at one level
                    %before the leaf level
                    else
                        %The current sub-cell will be subdivided 
                        %further, regardless of its control point signs
                        subcell_coords_all{lvl_d, occ_cell}((zcc_coords_cntr:(zcc_coords_cntr + 7)), 1:3) = sc_coords_orig;
                        zcc_coords_cntr = zcc_coords_cntr + 8;
                    end
                end %End sub_cell_8set      
            end %End lvl_d
        end %End check if lvl == (b + 1)  
        disp('------------------------------------------------------------');
    end %End occ_cell

    %Collect all of the sub-cell coordinates stored at the leaf level
    %of subcell_coords_all, for each occ_cell at level "lvl": these 
    %represent our reconstructed voxel (x, y, z) corner coordinates 
    %obtained from the reconstructed control points at level "lvl"
    disp(['Collecting all reconstructed voxels for each occupied cell at level ' num2str(lvl) ' and finding their centre coordinates ...']);
    reconstructed_vox_pos_corners = cat(1, subcell_coords_all{end, :});
    %Our input point cloud's voxel coordinates were considered to be 
    %the centres of the 1x1x1 voxels, so find the midpoint of each 
    %voxel in reconstructed_vox_pos_corners: these midpoints will 
    %represent our reconstructed voxel positions
    reconstructed_vox_pos = zeros((size(reconstructed_vox_pos_corners, 1)/8), 3);
    vox_cntr = 1;
    for vc = 1:8:(size(reconstructed_vox_pos_corners, 1) - 7)
        %Get the current set of 8 voxel corner coordinates
        vc_coords = reconstructed_vox_pos_corners((vc:(vc + 7)), :);
        %Find the mean of each of the 8 corner coordinates (x, y, and z
        %separately): these mean values represent the centre (x, y, z)
        %location of the current voxel
        reconstructed_vox_pos(vox_cntr, :) = mean(vc_coords, 1);
        vox_cntr = vox_cntr + 1;
    end

    %Order the reconstructed voxels according to their Morton codes, so
    %that they are in the same order as the input point cloud at the
    %encoder
    disp('Reordering reconstructed voxels according to Morton codes ...');
    %Get Morton codes for the reconstructed voxel x, y, z coordinates
    mortonCodes = xyzToMorton(reconstructed_vox_pos, lvl);   %"lvl" bits for each Morton code
    disp('Morton codes computed');
    %Sort the Morton codes obtained above, in ascending order
    [~, I_vox] = sort(mortonCodes);
    disp('Morton codes sorted');
    %Sort the voxel x, y, z locations in the same order as the sorted 
    %Morton codes
    reconstructed_vox_pos = reconstructed_vox_pos(I_vox, 1:3);
    disp('Reconstructed voxels sorted');
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
    %If we are not at the leaf level ...
    if lvl ~= (b + 1)
        title({'Voxel Reconstruction Using Only Reconstructed Control Points', ['at Octree Level ' num2str(lvl) ' and Interpolating at Higher Levels']});
    %If we are at the leaf level ...
    else
        title({'Voxel Reconstruction Using Only Reconstructed Control Points', ['at Octree Level ' num2str(lvl) ' (Leaf Level) - No Interpolation Required']});
    end
    %Save the above reconstruction as a MATLAB figure and as a PDF image in
    %our network directory (NB: The '-bestfit' option maximizes the size of 
    %the figure to fill the page, but preserves the aspect ratio of the 
    %figure. The figure might not fill the entire page. This option leaves 
    %a minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_cplvl' num2str(lvl)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_cplvl' num2str(lvl)], '-dpdf');
    disp('Saving reconstructed voxels figure ...');
    disp('------------------------------------------------------------');

    %For debugging purposes ...
    if lvl == (b + 1)
        %For debugging purposes, if we are at the leaf level, check if 
        %there are any voxels in the original voxelized point cloud 
        %that have not been reconstructed (i.e., are not present in 
        %reconstructed_vox_pos), and if so then plot these 
        test_vox_diffs = setdiff(ptcloud(:, 1:3), reconstructed_vox_pos, 'rows');
        disp(['Number of missing voxels in leaf level reconstruction: ' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' (' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']);
        if ~isempty(test_vox_diffs)
            figure;
            scatter3(test_vox_diffs(:, 1), test_vox_diffs(:, 2), test_vox_diffs(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
            axis equal; axis off;
            title({'Voxels that were Not Reconstructed at Leaf Level', ['(' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' = ' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']});
            savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\missing_voxels_leaf_level']);
            print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\missing_voxels_leaf_level'], '-dpdf');
            disp('Saving missing voxels figure ...');
            disp('------------------------------------------------------------');
        end
        %For debugging purposes, if we are at the leaf level, also 
        %check if any voxels are present in reconstructed_vox_pos that 
        %were NOT present in the original voxelized point cloud, and if
        %so then plot these
        test_vox_diffs2 = setdiff(reconstructed_vox_pos, ptcloud(:, 1:3), 'rows');
        disp(['Number of incorrectly reconstructed voxels at leaf level: ' num2str(size(test_vox_diffs2, 1))]);
        if ~isempty(test_vox_diffs2)
            figure;
            scatter3(test_vox_diffs2(:, 1), test_vox_diffs2(:, 2), test_vox_diffs2(:, 3), 5, 'filled', 'MarkerFaceColor', 'r');
            axis equal; axis off;
            title(['Incorrectly Reconstructed Voxels at Leaf Level: ' num2str(size(test_vox_diffs2, 1))]);
            savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_leaf_level']);
            print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_leaf_level'], '-dpdf');
            disp('Saving incorrect voxels figure ...');
            disp('------------------------------------------------------------');
        end
    end %End check if lvl == (b + 1) for debugging purposes    
end %End "lvl"
