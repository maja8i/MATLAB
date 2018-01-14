%For each leaf cell (Bezier Volume), find the dominant direction of the
%surface passing through that BV, based on the values of the control points
%at the corners of the BV, and ray-cast along this dominant direction from
%the midpoints of the 1x1 pixel positions on one of the faces (on the cube)
%that is orthogonal to that dominant direction.
function reconstructed_vox_pos = voxel_reconstruction_pruning_raycasting(debug_flag, pp_first_nonempty, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, q_stepsize)

%Matrix to store the reconstructed voxel positions (midpoints): initialize
%to some large value (can shrink later)
reconstructed_vox_pos = zeros(24000000, 3);
%Counter for reconstructed voxels
vox_cntr = 1;

%If there are any voxels that have already been reconstructed (because no
%pruning was done further up the octree for these voxels), get their corner
%coordinates and find their midpoints
if ~isempty(corner_coords_decoder{b + 1})
    if debug_flag == 1
        disp('Getting corner coordinates for already-reconstructed occupied voxels ...');
        disp(' ');    
    end
    reconstructed_vox_pos_corners = corner_coords_decoder{b + 1};
    for vc = 1:8:(size(reconstructed_vox_pos_corners, 1) - 7)
        %Get the current set of 8 voxel corner coordinates
        vc_coords = reconstructed_vox_pos_corners((vc:(vc + 7)), :);
        %Find the mean of each of the 8 corner coordinates (x, y, and z
        %separately): these mean values represent the centre (x, y, z) 
        %location of the current voxel
        reconstructed_vox_pos(vox_cntr, :) = mean(vc_coords, 1);
        vox_cntr = vox_cntr + 1;
    end
end

%Cell array that will store the indices of only the occupied cells from
%post_pruning_array that are leaves (i.e., our Bezier Volumes)
pp_leaves_only = cell(b, 1);
for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    %Find the location(s) of the leaf cell(s)
    leaves = find(post_pruning_array{lvl} == 1);
    %Store the indices of the leaf cells at this level
    pp_leaves_only{lvl} = leaves;
end

%For each octree level at which there are leaf cells, except the voxel
%level ...
for lvl = pp_first_nonempty:size(post_pruning_array, 1)
    %Find the side length of each leaf cell (BV) at the current level
    BV_side_length = (2^b)/(2^(lvl - 1));
    %For each leaf cell at this level ...
    for occ_cell = pp_leaves_only{lvl}'
        if debug_flag == 1
            disp(' ');
            disp(['Reconstructing voxels for leaf cell (occ_cell) ' num2str(occ_cell) ' at level ' num2str(lvl) ':']);
        end
        %Get the 8 corner coordinates of the current leaf cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
        %Get the control points for all 8 corners of the current leaf cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Find the dominant surface direction in this leaf cell, based on 
        %the values of its control points: this is the direction in which 
        %the partial derivative of the trilinear interpolation function has
        %the highest absolute value
        pd = zeros(3, 1);   %Partial derivatives vector
        pd(1) = current_ctrlpts(2) - current_ctrlpts(1) + current_ctrlpts(3) - current_ctrlpts(4) + current_ctrlpts(6) - current_ctrlpts(5) + current_ctrlpts(7) - current_ctrlpts(8);    %Partial derivative in x direction
        pd(2) = current_ctrlpts(4) - current_ctrlpts(1) + current_ctrlpts(3) - current_ctrlpts(2) + current_ctrlpts(7) - current_ctrlpts(6) + current_ctrlpts(8) - current_ctrlpts(5);    %Partial derivative in y direction
        pd(3) = current_ctrlpts(5) - current_ctrlpts(1) + current_ctrlpts(6) - current_ctrlpts(2) + current_ctrlpts(7) - current_ctrlpts(3) + current_ctrlpts(8) - current_ctrlpts(4);    %Partial derivative in z direction
        [max_pd, dominant_dir] = max(abs(pd));
        switch dominant_dir
            %Dominant direction is x
            case 1
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the x-axis 
                %(i.e., is on the face at x = 0, or on the face at x = 1)
                min_y = min(current_corner_coords(:, 2));
                max_y = max(current_corner_coords(:, 2));
                y_coords_norm = ((min_y:max_y) - min_y)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula
                y_coords = y_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                min_z = min(current_corner_coords(:, 3));
                max_z = max(current_corner_coords(:, 3));
                z_coords_norm = ((min_z:max_z) - min_z)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula
                z_coords = z_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                %Solve for x coordinates, to find any occupied voxels in
                %the current BV ...
                x_coords = zeros(length(y_coords)*length(z_coords), 1);  
                divisor_signs = [1;
                                -1;
                                -1;
                                 1;
                                 1;
                                -1;
                                -1;
                                 1];     
                %Solve the trilinear interpolation function f(x, y, z) = 0, 
                %for x, for each of the (y_coords, z_coords) pairs above
                vox_ind = 1;
                for i = 1:length(y_coords)     
                    for j = 1:length(z_coords)
                        divisor_yz_coords = [(1 - y_coords(i))*(1 - z_coords(j));
                                             (1 - y_coords(i))*(1 - z_coords(j));
                                                     y_coords(i)*(1-z_coords(j));
                                                     y_coords(i)*(1-z_coords(j));
                                                     (1-y_coords(i))*z_coords(j);
                                                     (1-y_coords(i))*z_coords(j);
                                                         y_coords(i)*z_coords(j);
                                                         y_coords(i)*z_coords(j)];

                        x_coords_divisor_temp = current_ctrlpts.*divisor_yz_coords.*divisor_signs;
                        x_coords_divisor = sum(x_coords_divisor_temp);
                        x_coords(vox_ind) = (current_ctrlpts(1)*(1-y_coords(i))*(1-z_coords(j)) + current_ctrlpts(4)*y_coords(i)*(1-z_coords(j)) + current_ctrlpts(5)*(1-y_coords(i))*z_coords(j) + current_ctrlpts(8)*y_coords(i)*z_coords(j))/x_coords_divisor;
                        vox_ind = vox_ind + 1;
                    end
                end  
                %Check which of the x_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                inds = find((x_coords >= 1/(2*BV_side_length)) & (x_coords <= (1 - 1/(2*BV_side_length))));
                [inds_r, inds_c] = ind2sub([length(y_coords) length(z_coords)], inds);
                occ_vox_x_coords = x_coords(inds);
                occ_vox_y_coords = y_coords(inds_r);
                occ_vox_z_coords = z_coords(inds_c);
                %Transform occupied voxels' x, y, z coordinates back from
                %[0, 1] to their original ranges
                occ_vox_y_coords = ((occ_vox_y_coords - (1/(2*BV_side_length))).*BV_side_length + min_y + 0.5)';
                occ_vox_z_coords = ((occ_vox_z_coords - (1/(2*BV_side_length))).*BV_side_length + min_z + 0.5)';
                min_x_midpoint = min(current_corner_coords(:, 1)) + 0.5;
                occ_vox_x_coords = (occ_vox_x_coords.*BV_side_length) + min_x_midpoint;
                %Round the x coordinates of the occupied voxels, to make
                %them integers
                occ_vox_x_coords = floor(occ_vox_x_coords);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];  
                vox_cntr = vox_cntr + length(occ_vox_x_coords);
            %Dominant direction is y
            case 2
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the y-axis 
                %(i.e., is on the face at y = 0, or on the face at y = 1)
                min_x = min(current_corner_coords(:, 1));
                max_x = max(current_corner_coords(:, 1));
                x_coords_norm = ((min_x:max_x) - min_x)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula
                x_coords = x_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                min_z = min(current_corner_coords(:, 3));
                max_z = max(current_corner_coords(:, 3));
                z_coords_norm = ((min_z:max_z) - min_z)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula
                z_coords = z_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                %Solve for y coordinates, to find any occupied voxels in
                %the current BV ...
                y_coords = zeros(length(x_coords)*length(z_coords), 1);
                divisor_signs = [1;
                                 1;
                                -1;
                                -1;
                                 1;
                                 1;
                                -1;
                                -1];
                %Solve the trilinear interpolation function f(x, y, z) = 0, 
                %for y, for each of the (x_coords, z_coords) pairs above 
                vox_ind = 1;
                for i = 1:length(x_coords)
                    for j = 1:length(z_coords)
                        divisor_xz_coords = [(1 - x_coords(i))*(1 - z_coords(j));
                                                   x_coords(i)*(1 - z_coords(j));
                                                   x_coords(i)*(1 - z_coords(j));
                                                 (1-x_coords(i))*(1-z_coords(j));
                                                     (1-x_coords(i))*z_coords(j);
                                                         x_coords(i)*z_coords(j);
                                                         x_coords(i)*z_coords(j);
                                                     (1-x_coords(i))*z_coords(j)];

                        y_coords_divisor_temp = current_ctrlpts.*divisor_xz_coords.*divisor_signs;
                        y_coords_divisor = sum(y_coords_divisor_temp);
                        y_coords(vox_ind) = (current_ctrlpts(1)*(1-x_coords(i))*(1-z_coords(j)) + current_ctrlpts(2)*x_coords(i)*(1-z_coords(j)) + current_ctrlpts(5)*(1-x_coords(i))*z_coords(j) + current_ctrlpts(6)*x_coords(i)*z_coords(j))/y_coords_divisor;
                        vox_ind = vox_ind + 1;
                    end
                end
                %Check which of the y_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                inds = find((y_coords >= 1/(2*BV_side_length)) & (y_coords <= (1 - 1/(2*BV_side_length))));
                [inds_r, inds_c] = ind2sub([length(x_coords) length(z_coords)], inds);
                occ_vox_y_coords = y_coords(inds);
                occ_vox_x_coords = x_coords(inds_r);
                occ_vox_z_coords = z_coords(inds_c);
                %Transform occupied voxels' x, y, z coordinates back from
                %[0, 1] to their original ranges
                occ_vox_x_coords = ((occ_vox_x_coords - (1/(2*BV_side_length))).*BV_side_length + min_x + 0.5)';
                occ_vox_z_coords = ((occ_vox_z_coords - (1/(2*BV_side_length))).*BV_side_length + min_z + 0.5)';
                min_y_midpoint = min(current_corner_coords(:, 2)) + 0.5;
                occ_vox_y_coords = (occ_vox_y_coords.*BV_side_length) + min_y_midpoint;
                %Round the y coordinates of the occupied voxels, to make
                %them integers
                occ_vox_y_coords = floor(occ_vox_y_coords);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];   
                vox_cntr = vox_cntr + length(occ_vox_x_coords);
            %Dominant direction is z    
            case 3
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the z-axis 
                %(i.e., is on the face at z = 0, or on the face at z = 1)
                min_x = min(current_corner_coords(:, 1));
                max_x = max(current_corner_coords(:, 1));
                x_coords_norm = ((min_x:max_x) - min_x)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula
                x_coords = x_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                min_y = min(current_corner_coords(:, 2));
                max_y = max(current_corner_coords(:, 2));
                y_coords_norm = ((min_y:max_y) - min_y)/BV_side_length;  %Coordinates have to be in [0, 1] for trilinear interpolation formula 
                y_coords = y_coords_norm(1:(end - 1)) + 1/(2*BV_side_length);
                %Solve for z coordinates, to find any occupied voxels in
                %the current BV ...
                z_coords = zeros(length(x_coords)*length(y_coords), 1);
                divisor_signs = [1;
                                 1;
                                 1;
                                 1;
                                -1;
                                -1;
                                -1;
                                -1];
                %Solve the trilinear interpolation function f(x, y, z) = 0, 
                %for z, for each of the (x_coords, y_coords) pairs above
                vox_ind = 1;
                for i = 1:length(x_coords)
                    for j = 1:length(y_coords)
                        divisor_xy_coords = [(1 - x_coords(i))*(1 - y_coords(j));
                                                   x_coords(i)*(1 - y_coords(j));
                                                         x_coords(i)*y_coords(j);
                                                     (1-x_coords(i))*y_coords(j);
                                                 (1-x_coords(i))*(1-y_coords(j));
                                                     x_coords(i)*(1-y_coords(j));
                                                         x_coords(i)*y_coords(j);
                                                     (1-x_coords(i))*y_coords(j)];

                        z_coords_divisor_temp = current_ctrlpts.*divisor_xy_coords.*divisor_signs;
                        z_coords_divisor = sum(z_coords_divisor_temp);
                        z_coords(vox_ind) = (current_ctrlpts(1)*(1-x_coords(i))*(1-y_coords(j)) + current_ctrlpts(2)*x_coords(i)*(1-y_coords(j)) + current_ctrlpts(3)*x_coords(i)*y_coords(j) + current_ctrlpts(4)*(1-x_coords(i))*y_coords(j))/z_coords_divisor;
                        vox_ind = vox_ind + 1;
                    end
                end
                %Check which of the z_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                inds = find((z_coords >= 1/(2*BV_side_length)) & (z_coords <= (1 - 1/(2*BV_side_length))));
                [inds_r, inds_c] = ind2sub([length(x_coords) length(y_coords)], inds);
                occ_vox_z_coords = z_coords(inds);
                occ_vox_x_coords = x_coords(inds_r);
                occ_vox_y_coords = y_coords(inds_c);
                %Transform occupied voxels' x, y, z coordinates back from
                %[0, 1] to their original ranges
                occ_vox_x_coords = ((occ_vox_x_coords - (1/(2*BV_side_length))).*BV_side_length + min_x + 0.5)';
                occ_vox_y_coords = ((occ_vox_y_coords - (1/(2*BV_side_length))).*BV_side_length + min_y + 0.5)';
                min_z_midpoint = min(current_corner_coords(:, 3)) + 0.5;
                occ_vox_z_coords = (occ_vox_z_coords.*BV_side_length) + min_z_midpoint;
                %Round the z coordinates of the occupied voxels, to make
                %them integers
                occ_vox_z_coords = floor(occ_vox_z_coords);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];    
                vox_cntr = vox_cntr + length(occ_vox_x_coords);
        end %End switch max_pd
    end %End occ_cell
end %End lvl

%If there are any surplus rows in reconstructed_vox_pos (because the
%initial matrix size was made too large), remove them now
reconstructed_vox_pos((vox_cntr:end), :) = [];




