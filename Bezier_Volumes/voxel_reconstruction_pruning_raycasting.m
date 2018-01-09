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
        pd(1) = (current_ctrlpts(2) - current_ctrlpts(1)) + (current_ctrlpts(3) - current_ctrlpts(4)) + (current_ctrlpts(6) - current_ctrlpts(5)) + (current_ctrlpts(7) - current_ctrlpts(8));    %Partial derivative in x direction
        pd(2) = (current_ctrlpts(4) - current_ctrlpts(1)) + (current_ctrlpts(3) - current_ctrlpts(2)) + (current_ctrlpts(7) - current_ctrlpts(6)) + (current_ctrlpts(8) - current_ctrlpts(5));    %Partial derivative in y direction
        pd(3) = (current_ctrlpts(5) - current_ctrlpts(1)) + (current_ctrlpts(6) - current_ctrlpts(2)) + (current_ctrlpts(7) - current_ctrlpts(3)) + (current_ctrlpts(8) - current_ctrlpts(4));    %Partial derivative in z direction
        [max_pd, dominant_dir] = max(abs(pd));
        switch max_pd
            %Dominant direction is x
            case 1
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the x-axis 
                %(i.e., is on the face at x = 0, or on the face at x = 1)
                y_coords = ;    %Coordinates have to be in the range [0, 1] for the trilinear interpolation formula to work
                z_coords = ;
                x_coords = zeros(length(y_coords), 1);
                
                divisor_signs = [1;
                                -1;
                                -1;
                                 1;
                                 1;
                                -1;
                                -1;
                                 1];
                             
                for i = 1:length(y_coords)
                    %Solve the trilinear interpolation function f(x, y, z) = 0, 
                    %for x, for each of the (y_coords, z_coords) pairs above     
                    divisor_yz_coords = [(1 - y_coords(i))*(1 - z_coords(i));
                                         (1 - y_coords(i))*(1 - z_coords(i));
                                                 y_coords(i)*(1-z_coords(i));
                                                 y_coords(i)*(1-z_coords(i));
                                                 (1-y_coords(i))*z_coords(i);
                                                 (1-y_coords(i))*z_coords(i);
                                                     y_coords(i)*z_coords(i);
                                                     y_coords(i)*z_coords(i)];

                    x_coords_divisor_temp = current_ctrlpts.*divisor_yz_coords.*divisor_signs;
                    x_coords_divisor = sum(x_coords_divisor_temp);

                    x_coords(i) = (current_ctrlpts(1)*(1-y_coords(i))*(1-z_coords(i)) + current_ctrlpts(4)*y_coords(i)*(1-z_coords(i)) + current_ctrlpts(5)*(1-y_coords(i))*z_coords(i) + current_ctrlpts(8)*y_coords(i)*z_coords(i))/x_coords_divisor;
                end  
                %Check which of the x_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                [occ_vox_x_coords, inds] = find(((x_coords >= 0) && (x_coords <= 1)));
                occ_vox_y_coords = y_coords(inds);
                occ_vox_z_coords = z_coords(inds);
                %Quantize the x_coords of the occupied voxels found above
                occ_vox_x_coords = quantize_uniform_scalar(occ_vox_x_coords, q_stepsize);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];  
                vox_cntr = vox_cntr + size(occ_vox_x_coords, 1);
            %Dominant direction is y
            case 2
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the y-axis 
                %(i.e., is on the face at y = 0, or on the face at y = 1)
                x_coords = ;    %Coordinates have to be in the range [0, 1] for the trilinear interpolation formula to work
                z_coords = ;
                y_coords = zeros(length(x_coords), 1);
                
                divisor_signs = [1;
                                 1;
                                -1;
                                -1;
                                 1;
                                 1;
                                -1;
                                -1];
              
                for i = 1:length(x_coords)
                    %Solve the trilinear interpolation function f(x, y, z) = 0, 
                    %for y, for each of the (x_coords, z_coords) pairs above                    
                    divisor_xz_coords = [(1 - x_coords(i))*(1 - z_coords(i));
                                               x_coords(i)*(1 - z_coords(i));
                                               x_coords(i)*(1 - z_coords(i));
                                             (1-x_coords(i))*(1-z_coords(i));
                                                 (1-x_coords(i))*z_coords(i);
                                                     x_coords(i)*z_coords(i);
                                                     x_coords(i)*z_coords(i);
                                                 (1-x_coords(i))*z_coords(i)];

                    y_coords_divisor_temp = current_ctrlpts.*divisor_xz_coords.*divisor_signs;
                    y_coords_divisor = sum(y_coords_divisor_temp);

                    y_coords(i) = (current_ctrlpts(1)*(1-x_coords(i))*(1-z_coords(i)) + current_ctrlpts(2)*x_coords(i)*(1-z_coords(i)) + current_ctrlpts(5)*(1-x_coords(i))*z_coords(i) + current_ctrlpts(6)*x_coords(i)*z_coords(i))/y_coords_divisor;
                end
                %Check which of the y_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                [occ_vox_y_coords, inds] = find(((y_coords >= 0) && (y_coords <= 1)));
                occ_vox_x_coords = x_coords(inds);
                occ_vox_z_coords = z_coords(inds);
                %Quantize the y_coords of the occupied voxels found above
                occ_vox_y_coords = quantize_uniform_scalar(occ_vox_y_coords, q_stepsize);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];   
                vox_cntr = vox_cntr + size(occ_vox_y_coords, 1);
            %Dominant direction is z    
            case 3
                %Find all 1x1 pixels that can fit on one face of the 
                %current BV (occ_cell) that is orthogonal to the z-axis 
                %(i.e., is on the face at z = 0, or on the face at z = 1)
                x_coords = ;    %Coordinates have to be in the range [0, 1] for the trilinear interpolation formula to work
                y_coords = ;
                z_coords = zeros(length(x_coords), 1);
                
                divisor_signs = [1;
                                 1;
                                 1;
                                 1;
                                -1;
                                -1;
                                -1;
                                -1];
                
                for i = 1:length(x_coords)
                    %Solve the trilinear interpolation function f(x, y, z) = 0, 
                    %for z, for each of the (x_coords, y_coords) pairs above
                    divisor_xy_coords = [(1 - x_coords(i))*(1 - y_coords(i));
                                               x_coords(i)*(1 - y_coords(i));
                                                     x_coords(i)*y_coords(i);
                                                 (1-x_coords(i))*y_coords(i);
                                             (1-x_coords(i))*(1-y_coords(i));
                                                 x_coords(i)*(1-y_coords(i));
                                                     x_coords(i)*y_coords(i);
                                                 (1-x_coords(i))*y_coords(i)];

                    z_coords_divisor_temp = current_ctrlpts.*divisor_xy_coords.*divisor_signs;
                    z_coords_divisor = sum(z_coords_divisor_temp);

                    z_coords(i) = (current_ctrlpts(1)*(1-x_coords(i))*(1-y_coords(i)) + current_ctrlpts(2)*x_coords(i)*(1-y_coords(i)) + current_ctrlpts(3)*x_coords(i)*y_coords(i) + current_ctrlpts(4)*(1-x_coords(i))*y_coords(i))/z_coords_divisor;
                end
                %Check which of the z_coords are in the range [0, 1]: these
                %represent occupied voxels on the surface inside the
                %current BV block
                [occ_vox_z_coords, inds] = find(((z_coords >= 0) && (z_coords <= 1)));
                occ_vox_x_coords = x_coords(inds);
                occ_vox_y_coords = y_coords(inds);
                %Quantize the z_coords of the occupied voxels found above
                occ_vox_z_coords = quantize_uniform_scalar(occ_vox_z_coords, q_stepsize);  
                %Record the occupied voxels' (x, y, z) positions (i.e., of
                %their midpoints)
                reconstructed_vox_pos(vox_cntr:(vox_cntr + length(occ_vox_x_coords) - 1), :) = [occ_vox_x_coords occ_vox_y_coords occ_vox_z_coords];    
                vox_cntr = vox_cntr + size(occ_vox_z_coords, 1);
        end %End switch max_pd
    end %End occ_cell
    
    
end %End lvl

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



