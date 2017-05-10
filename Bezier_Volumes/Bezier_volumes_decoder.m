%Requires directory "phil": add this directory plus its sub-directories to
%the current MATLAB path.

%start_lvl = 1;
%max_lvl = 8;
%b = 10;
%rec_ctrlpts_start_lvl = reconstructed_control_points{start_lvl};
%OccupancyCode = myOT.OccupancyCode; %Should be from start_lvl only
%vis_levels_ot = 4; %No. of octree levels for which we want to visualize the octree cell subdivision

function [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name)

disp(' ');
disp('============================================================');
disp('                   DECODER RUNNING ...');
disp('============================================================');
disp(' ');

OccupancyCode = occupancy_codes_forDec;
rec_ctrlpts_start_lvl = rec_ctrlpts_forDec;
wavelet_coeffs = wavelet_coeffs_forDec;

%------------------------- Octree Reconstruction -------------------------%

disp('------------------- Octree Reconstruction ------------------');
disp(' ');

tic;

%Initialize a cell array to store the number of children for each occupied
%octree cell at each level
ChildCount = cell(b, 1);

%Initialize a cell array to store the pointer to the first child of each 
%occupied octree cell at each level
FirstChildPtr = cell(b, 1);

%Initialize a cell array to store the SpatialIndex for each occupied octree
%cell at each level
SpatialIndex = cell((b + 1), 1);
%Initialize SpatialIndex at the root to [0, 0, 0]
SpatialIndex{1} = [uint16(0), uint16(0), uint16(0)];

%Initialize a matrix to serve as the SpatialIndex starting point
%(dictionary)
SI_dict = zeros(8, 3);
SI_dict(1, :) = SpatialIndex{1};    %Origin (0, 0, 0)
SI_dict(2, :) = [0 0 1];    %(+z)
SI_dict(3, :) = [0 1 0];    %(+y)
SI_dict(4, :) = [0 1 1];    %(+y, +z)
SI_dict(5, :) = [1 0 0];    %(+x)
SI_dict(6, :) = [1 0 1];    %(+x, +z)
SI_dict(7, :) = [1 1 0];    %(+x, +y)
SI_dict(8, :) = [1 1 1];    %(+x, +y, +z)

%For each octree level ...
%for lvl = 1:b
for lvl = 1:(max_lvl - 1)
    disp(['Processing octree level ' num2str(lvl) ':']);
    %Counter for all children of all occupied nodes at this level
    total_child_cntr = 1;
    %For each occupied cell at this level ...
    for occ_cell = 1:numel(OccupancyCode{lvl})
        %Convert the OccupancyCode decimal value for this cell's children,
        %into its binary representation
        bin_vec = dec2bin(OccupancyCode{lvl}(occ_cell), 8);
        %The number of "1"s in bin_vec indicates the number of children
        %that this octree cell has
        ChildCount{lvl}(occ_cell) = numel(strfind(bin_vec, '1'));    
        %Find the locations in bin_vec where the vector contains '1's
        ones_inds = strfind(bin_vec, '1');
        %Compute SpatialIndex for each occupied child    
        for occ_child = 1:length(ChildCount{lvl}(occ_cell))
            %If we are currently at the root cell, we can obtain its
            %children's spatial indices just by indexing into SI_dict
            if lvl == 1
                SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint16(SI_dict(ones_inds, :));
            else
                %Find where the curent occupied cell's spatial index 
                %contains non-zero values: these are the positions to which
                %we will add offsets to obtain the children's spatial
                %indices, when necessary
                parent_SI_nonzero = uint16(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                %If parent_SI_nonzero is empty (i.e., the current occupied
                %cell's spatial index is (0, 0, 0)), then the spatial
                %indices of this cell's children can be obtained just by
                %indexing into SI_dict
                if isempty(parent_SI_nonzero)
                    SpatialIndex{lvl + 1}(1:ChildCount{lvl}(occ_cell), :) = uint16(SI_dict(ones_inds, :));
                else
                    %Add an offset to corresponding locations in SI_dict,
                    %in the columns determined by parent_SI_nonzero
                    SI_dict_locations = SI_dict(ones_inds, :);
                    SI_dict_locations(:, parent_SI_nonzero) = uint16(SI_dict_locations(:, parent_SI_nonzero)) + uint16(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                    SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(occ_cell) - 1), :) = uint16(SI_dict_locations);
                end
            end
        end %End occ_child  
        total_child_cntr = total_child_cntr + ChildCount{lvl}(occ_cell);
    end %End occ_cell
    %Make ChildCount{lvl} a column vector rather than a row vector
    ChildCount{lvl} = ChildCount{lvl}';
    disp('Finished computing ChildCount for each occupied cell at this level, and SpatialIndex for each occupied child of each occupied cell');
    %Compute the first child pointer for each occupied cell at the current
    %octree level
    lastChildPtr = cumsum(int32(ChildCount{lvl}));
    FirstChildPtr{lvl} = uint32(lastChildPtr - int32(ChildCount{lvl}) + 1);
    disp('Finished computing FirstChildPtr for each occupied cell at this level');
    disp('------------------------------------------------------------');
end %End lvl
OT_recon_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken for octree reconstruction: ' num2str(OT_recon_time) ' seconds']);
disp('************************************************************');

%--------------------- Computing Corner Coordinates ----------------------%

disp(' ');
disp('-------------- Computing Corner Coordinates ----------------');
disp(' ');

%Initialize a cell array to store the corner coordinates of each occupied
%cell at each level of the octree
corner_coords_decoder = cell((b + 1), 1);
%Initialize a cell array to store only the unique corner coordinates at
%each octree level
unique_coords_decoder = cell((b + 1), 1);
%Initialize a cell array of "pointers": there will be 8 pointers per
%occupied octree cell (1 per corner) at each octree level, which will point
%to the Bezier control point in control_points, associated with that
%corner. So, for each octree level "lvl", there will be 
%numel(OccupancyCode{lvl}) pointer arrays containining 8 elements each.  
ctrl_pts_pointers = cell((b + 1), 1);

%Figure out the corner coordinates for each corner of each occupied cell at
%each octree level, starting from start_lvl
tic;
%for lvl = start_lvl:(b + 1)
for lvl = start_lvl:max_lvl
%     %Counter to keep track of how many corner coordinates we have stored at
%     %the current octree level
%     corner_coords_cntr = 1;
%     %For each occupied cell at the current level
%     for occ_cell = 1:numel(OccupancyCode{lvl})
%         %Initialize a matrix to store the (x, y, z) coordinates of all the
%         %corners of the current octree cell
%         corners = zeros(8, 3);
%         %Find the (x, y, z) coordinates of the origin of this cell (found
%         %at the bottom left-hand corner farthest from the viewer)
%         corners(1, :) = double(SpatialIndex{lvl}(occ_cell, :))*(2^(b + 1 - lvl)) - [0.5 0.5 0.5];
%         %Find the (x, y, z) coordinates of the other 7 corners of this cell
%         corners(2, :) = corners(1, :) + [2^(b + 1 - lvl) 0 0];
%         corners(3, :) = corners(1, :) + [2^(b + 1 - lvl) 2^(b + 1 - lvl) 0];
%         corners(4, :) = corners(1, :) + [0 2^(b + 1 - lvl) 0];
%         corners(5, :) = corners(1, :) + [0 0 2^(b + 1 - lvl)];
%         corners(6, :) = corners(1, :) + [2^(b + 1 - lvl) 0 2^(b + 1 - lvl)];
%         corners(7, :) = corners(1, :) + [2^(b + 1 - lvl) 2^(b + 1 - lvl) 2^(b + 1 - lvl)];
%         corners(8, :) = corners(1, :) + [0 2^(b + 1 - lvl) 2^(b + 1 - lvl)];
%         %Store all the corner coordinates for the current cell in their
%         %corresponding locations inside corner_coords
%         corner_coords_decoder{lvl}(corner_coords_cntr:(corner_coords_cntr + 7), 1:3) = corners;   
%         corner_coords_cntr = corner_coords_cntr + 8;
%     end

    disp(['Processing octree level ' num2str(lvl) ':']); 
    disp('Computing corner coordinates for each occupied cell ...');
    disp('------------------------------------------------------------');
    
    %Find the (x, y, z) coordinates of the origin of each occupied octree
    %cell at the current level (origin is at the bottom left-hand corner 
    %farthest from the viewer)
    corners1 = double(SpatialIndex{lvl})*(2^(b + 1 - lvl)) - [0.5 0.5 0.5];
    %Replicate each row of corners1 7 times, so that we can directly add
    %this matrix to offsets_from_origin_rep (below)
    corners1_rep = corners1(repmat(1:size(corners1, 1), 7, 1), :);
    
    %Find the (x, y, z) coordinates of the other 7 corners of each occupied
    %octree cell at the current level
    offsets_from_origin = [[2^(b + 1 - lvl) 0 0]; [2^(b + 1 - lvl) 2^(b + 1 - lvl) 0]; [0 2^(b + 1 - lvl) 0]; [0 0 2^(b + 1 - lvl)]; [2^(b + 1 - lvl) 0 2^(b + 1 - lvl)]; [2^(b + 1 - lvl) 2^(b + 1 - lvl) 2^(b + 1 - lvl)]; [0 2^(b + 1 - lvl) 2^(b + 1 - lvl)]];
    offsets_from_origin_rep = repmat(offsets_from_origin, size(corners1, 1), 1);
    corners2_8 = corners1_rep + offsets_from_origin_rep;
    
    %Store all the corner coordinates for the current cell in their
    %corresponding locations inside corner_coords
    corner_coords_decoder{lvl} = zeros((size(SpatialIndex{lvl}, 1)*8), 3);
    corner_coords_decoder{lvl}(1:8:(size(SpatialIndex{lvl}, 1)*8 - 7), :) = corners1;   
    corners2_8_cntr = 1;
    for next_ind = 2:8:(size(SpatialIndex{lvl}, 1)*8 - 6)
        corner_coords_decoder{lvl}((next_ind:(next_ind + 6)), :) = corners2_8((corners2_8_cntr:(corners2_8_cntr + 6)), :);
        corners2_8_cntr = corners2_8_cntr + 7;
    end

    %Find all the unique corner coordinates at the current level of the
    %octree, and in ctrl_pts_pointers{lvl}, for each row that represents an
    %(x, y, z) coordinate from corner_coords{lvl}, store a pointer for this
    %coordinate into the unique_coords{lvl} array, to say what the
    %coordinate triplet in this row should be. These pointers will also act
    %as pointers to the reconstructed control points stored in 
    %reconstruction_decoder{lvl}, since reconstruction_decoder{lvl} will be 
    %the same size as unique_coords{lvl}, as each unique corner coordinate 
    %triplet will have one corresponding Bezier control point (signed 
    %distance value) associated with it.
    [unique_coords_decoder{lvl}, ~, ctrl_pts_pointers{lvl}] = unique(corner_coords_decoder{lvl}, 'rows', 'stable');
end
cornercoords_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all corner coordinates: ' num2str(cornercoords_time) ' seconds']);
disp('************************************************************');

%---------------- Signal (Control Point) Reconstruction ------------------%

disp(' ');
disp('-------------- Control Point Reconstruction ----------------');
disp(' ');

%Initialize a cell array to store the reconstructed signal (control point)
%at each vertex of each octree cell at every level from start_lvl to the 
%leaves. Note that we can't just use the reconstructed_control_points from 
%the encoder, as we won't be transmitting those (except for the 
%reconstructed control points at start_lvl)
reconstruction_decoder = cell((b + 1), 1);

%Dequantize the reconstructed control points for start_lvl, which were 
%received from the encoder, then place them in the corresponding cell of 
%reconstruction_decoder. Note that when q_stepsize = 1, the dequantized
%values will just be equal to the quantized values.
reconstruction_decoder{start_lvl} = dequantize_uniform_scalar(rec_ctrlpts_start_lvl, q_stepsize);

%Start from the start_lvl and work up to the leaves ...
tic;
%for lvl = start_lvl:b
for lvl = start_lvl:(max_lvl - 1)
    disp(['Reconstructing control points at octree level ' num2str(lvl + 1) ' ...']);
    disp('------------------------------------------------------------');
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied cell at the current octree level ...
    for occ_cell = 1:numel(OccupancyCode{lvl})
        %Extract the cell's 8 corner coordinates. This cell will represent
        %our parent cell at the current level, since we will consider its
        %children for the control point (signal) reconstruction.
        parent_corner_coords = corner_coords_decoder{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);
        %Get the pointer to the first child cell of the current parent cell
        child_ptr = FirstChildPtr{lvl}(occ_cell);
        %For each child of the current cell ...
        for child = 1:ChildCount{lvl}(occ_cell)
            %Get the 8 corner coordinates of the current child cell. NOTE:
            %Below, "child" is converted to type double, because it is
            %uint8 by default, and if the result of any of the 
            %multiplications is > 255, the answer will be truncated and the
            %final result will be incorrect.
            child_corner_coords = corner_coords_decoder{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
            %For each corner of the current child cell ...
            for cnr = 1:8
                %Check if the current corner already has a reconstructed 
                %control point associated with it (since some corners will  
                %be shared amongst different octree cells and we want to 
                %make sure that we process only the UNIQUE corner vertices 
                %at each octree level) 
                possible_inds = ctrl_pts_pointers{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
                if ~isempty(reconstruction_decoder{lvl + 1})
                    if length(reconstruction_decoder{lvl + 1}) >= possible_inds(cnr) 
                        if ~isempty(reconstruction_decoder{lvl + 1}(possible_inds(cnr)))
                            continue;          
                        end
                    end
                end 
                %Flag to indicate if current corner is on a parent edge 
                %(0 => no; 1 => yes)
                on_p_edge = 0;
                %Flag to indicate if current corner is on a parent face 
                %(0 => no; 1 => yes)
                on_p_face = 0;
                %For the current corner, check if it is on a parent edge 
                %(i.e., if it shares at least 2 same coordinates (out of x,
                %y, or z) with two of the parent vertices), or on a parent 
                %face (i.e., only 1 of its coordinates (either x, or y, or 
                %z) is the same as four parents' coordinates (must have the 
                %same coordinate in common in this case: either x, or y, or 
                %z))
                parent_row_inds = [];
                [x_same, ~] = find(child_corner_coords(cnr, 1) == parent_corner_coords(:, 1));
                [y_same, ~] = find(child_corner_coords(cnr, 2) == parent_corner_coords(:, 2));
                [z_same, ~] = find(child_corner_coords(cnr, 3) == parent_corner_coords(:, 3));
                temp_cat = [x_same; y_same; z_same];
                if ((isempty(x_same) + isempty(y_same) + isempty(z_same)) == 2)
                    %Corner is on a parent face
                    on_p_face = 1;
                    %Get the row indices of the parent vertices on this face
                    if ~isempty(x_same)
                        parent_row_inds = x_same;
                    elseif ~isempty(y_same)
                        parent_row_inds = y_same;
                    elseif ~isempty(z_same)
                        parent_row_inds = z_same;
                    end
                else
                    %Check if any vertex indices appear twice in temp_cat
                    temp_cntr = 1;
                    for v = 1:length(temp_cat)
                        if length(find(temp_cat == temp_cat(v))) == 2
                            parent_row_inds(temp_cntr) = temp_cat(v);
                            temp_cntr = temp_cntr + 1;
                        end
                    end
                    if (length(unique(parent_row_inds)) == 2)
                        %Corner is on a parent edge
                        on_p_edge = 1;
                        %Keep only the unique indices in parent_row_inds
                        parent_row_inds = unique(parent_row_inds);
                    end      
                end
                
                %If this corner's coordinates are the same as one of the
                %corner coordinates of the parent
                if sum(ismember(parent_corner_coords, child_corner_coords(cnr, :), 'rows') > 0)
                    %Get the index of the current parent corner
                    %parent_row_index = find(ismember(corner_coords_decoder{lvl}, child_corner_coords(cnr, :), 'rows') > 0);
                    parent_row_index = find(sum(abs(child_corner_coords(cnr, :) - corner_coords_decoder{lvl}), 2) == 0, 1); %Faster than "ismember"
                    %Find the control point index for the parent corner.
                    %Although the length of parent_row_index may sometimes
                    %be > 1, and so more than one control point index may
                    %be found below, in this case the result should still 
                    %be the same index, just repeated. So extract only the 
                    %unique control point index found below (should be just
                    %one).
                    %parent_ctrlpt_index = unique(ctrl_pts_pointers{lvl}(parent_row_index));
                    parent_ctrlpt_index = ctrl_pts_pointers{lvl}(parent_row_index);
                    %Do nothing (the signal on this vertex is a low-pass
                    %coefficient and has already been reconstructed),
                    %except transfer the reconstructed control point over
                    %from the parent corner at the previous octree level 
                    %(we want the reconstruction_decoder values at each 
                    %level to correspond to the same locations in 
                    %unique_coords)
                    reconstruction_decoder{lvl + 1}(possible_inds(cnr)) = reconstruction_decoder{lvl}(parent_ctrlpt_index);
                    continue;
                %If this corner vertex lies on a parent edge
                elseif on_p_edge == 1
                    %Get the reconstructed Bezier control points stored at 
                    %the corner vertices defined by parent_row_inds
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt1 = reconstruction_decoder{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstruction_decoder{lvl}(ctrlpt2_ptr);
                    %Average the signal (Bezier control points) on the 2
                    %vertices of the parent edge
                    avg_signal = (ctrlpt1 + ctrlpt2)/2;
                %If this corner vertex lies on a parent face 
                elseif on_p_face == 1
                    %Get the reconstructed Bezier control points stored at 
                    %the corner vertices defined by parent_row_inds
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt3_ptr = all_ctrlpts_ptrs(parent_row_inds(3));
                    ctrlpt4_ptr = all_ctrlpts_ptrs(parent_row_inds(4));
                    ctrlpt1 = reconstruction_decoder{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstruction_decoder{lvl}(ctrlpt2_ptr);
                    ctrlpt3 = reconstruction_decoder{lvl}(ctrlpt3_ptr);
                    ctrlpt4 = reconstruction_decoder{lvl}(ctrlpt4_ptr);
                    %Average the signal (Bezier control points) on the 4 
                    %vertices of the parent face
                    avg_signal = (ctrlpt1 + ctrlpt2 + ctrlpt3 + ctrlpt4)/4;  
                %If this corner vertex lies somewhere in the centre of the
                %parent's block (i.e., neither on a parent's edge or on a
                %parent's face)
                else                
                    %Get the reconstructed Bezier control points stored at 
                    %all of the corner vertices (8) of the parent cell
                    ctrlpts_pointers = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpts = reconstruction_decoder{lvl}(ctrlpts_pointers);
                    %Average the signal (Bezier control points) on the 8
                    %vertices of the parent block
                    avg_signal = mean(ctrlpts);
                end
                %Get the dequantized wavelet coefficient for the current
                %corner of the current octree cell (the child cell)
                child_wavelet_coeff = dequantize_uniform_scalar(wavelet_coeffs{lvl + 1}(possible_inds(cnr)), q_stepsize);
                %Add the dequantized wavelet coefficient to avg_signal, to 
                %obtain the reconstructed signal (control point) at the 
                %current child vertex 
                reconstruction_decoder{lvl + 1}(possible_inds(cnr)) = child_wavelet_coeff + avg_signal;
            end %End corners 
        end %End children
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
    end %End occupied cells at current level
    %Arrange all the reconstructed control points produced for the child
    %octree level, into a column vector
    reconstruction_decoder{lvl + 1} = (reconstruction_decoder{lvl + 1})';
end %End "lvl"
ctrlpt_recon_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to reconstruct all control points: ' num2str(ctrlpt_recon_time) ' seconds']);
disp('************************************************************');

% %--------------- Reconstructed Bezier Volume Visualization ---------------%
% 
% disp(' ');
% disp('-------- Reconstructed Bezier Volume Visualization ---------');
% disp(' ');
% 
% %For each octree level for which we have reconstructed control points, plot
% %a sphere of possible nearest voxel points for each unique corner
% %coordinate (the possible nearest voxel points would be on the surface of
% %this sphere)
% %for lvl = 1:size(reconstruction_decoder, 1)
% for lvl = 1:6
%     disp(['Plotting nearest voxel spheres at octree level ' num2str(lvl) ' ...']);
%     disp('------------------------------------------------------------');
%     figure;
%     %For each unique corner at this octree level ...
%     for un_cnr = 1:size(unique_coords_decoder{lvl}, 1)
%         %Decoder does not know the exact coordinate values of the nearest
%         %(occupied) voxel point to this corner, but the absolute value of 
%         %the signed distance value corresponding to this corner represents 
%         %the radius of the sphere that represents all the possible nearest
%         %voxel points that could correspond to this corner (on the surface
%         %of this sphere)
%         h(1) = scatter3(unique_coords{lvl}(un_cnr, 1), unique_coords{lvl}(un_cnr, 2), unique_coords{lvl}(un_cnr, 3), 40, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); %Plot the current unique corner point
%         hold on;
%         %Draw a sphere centred at the current corner coordinate, with a 
%         %radius equal to the absolute value of the reconstructed control 
%         %point associated with this corner
%         [xs, ys, zs] = sphere;
%         surf(xs, ys, zs);
%         hold on;
%         h(2) = surf((abs(reconstruction_decoder{lvl}(un_cnr))*xs + unique_coords{lvl}(un_cnr, 1)), (abs(reconstruction_decoder{lvl}(un_cnr))*ys + unique_coords{lvl}(un_cnr, 2)), (abs(reconstruction_decoder{lvl}(un_cnr))*zs + unique_coords{lvl}(un_cnr, 3)));             
%         set(h(2), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceLighting', 'phong');  %Make sphere semi-transparent
%         axis equal;
%         view(-58, 39);
%     end
%     legend(h(1), 'Unique corner coordinate', 'Location', 'best');
%     title({'Spheres of Possible Nearest Voxel Coordinates', ['for Each Unique Corner Coordinate at Octree Level ' num2str(lvl)]});
% end

%------ Point Cloud Voxel Reconstruction at Different Octree Levels ------%

disp(' ');
disp('------------------ Voxel Reconstruction --------------------');
disp(' ');

%For a chosen (set of) octree level(s), reconstruct the voxel positions 
%that define the surface (shape) of the input 3D point cloud, by using only 
%the reconstructed control points at the chosen octree level(s) and 
%interpolating between them for the higher octree levels, to figure out 
%where all of the zero crossings are (i.e., through which octree cells the
%surface passes), down to the leaf level. Finally, at the leaf level, the 
%midpoints of the voxels that contain the zero crossings will be our 
%reconstructed voxel positions.
%for lvl = [4, 5, 6, 7, 8]
for lvl = max_lvl
    disp(['Reconstructing voxels for level ' num2str(lvl) ' control points ...']);
    disp(' ');
    subcell_coords_all = cell((b + 1), size(SpatialIndex{lvl}, 1));    %Cell array to store the coordinates of octree cells at all levels from (lvl + 1) to the leaves, which contain zero crossings
    %For each occupied cell at the current octree level, lvl ...
    for occ_cell = 1:size(SpatialIndex{lvl}, 1) %We can only use this for levels for which we know how many occupied cells there are
        disp(['Reconstructing voxels for occ_cell ' num2str(occ_cell) '/' num2str(size(SpatialIndex{lvl}, 1)) ':']);
        %Get the control points for all 8 corners of this cell
        current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %Get the corner coordinates for each of the 8 corners of this cell
        current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
        %Check if all 8 control points for this cell have the same sign 
        if (sum(sign(current_ctrlpts)) == length(current_ctrlpts))||(sum(sign(current_ctrlpts)) == -length(current_ctrlpts))
            %If yes, this indicates that the surface of the 3D object does not 
            %pass through this octree block (i.e., there is no zero crossing 
            %here), so this cell is not a candidate for further subdivision; do
            %nothing further 
            disp('No zero crossings in this cell');
            disp('------------------------------------------------------------');
            continue;  
        else 
            %If no (i.e., if there is a sign change within the cell), this 
            %indicates that the surface of the 3D object cuts through this cell
            %(i.e., there is a zero crossing in this cell), so this cell will 
            %need to be subdivided further if it's not at the leaf level 
            %(that is, we will need to look at 
            %the interpolated control points of all of this cell's children, 
            %grandchildren and other descendants, to check for zero crossings 
            %and thereby figure out which cells are occupied at each octree 
            %level, down to the leaves)
            if lvl == (b + 1)
                %Just get the current occupied cell's (voxel's) corner 
                %coordinates 
                subcell_coords_all{lvl, occ_cell} = current_corner_coords;
            else %lvl ~= max_lvl
                %For each descendant of the current occ_cell, right down to
                %the leaves ...
                for lvl_d = (lvl + 1):(b + 1)      
                    disp(['--Processing sub-cells at level ' num2str(lvl_d)]);
                    %Find the unique minimum, unique maximum, and midpoint
                    %coordinates for the x, y, and z dimensions of the sub-cells at 
                    %this level
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
                        %There will be one unq_min, one unq_max, and one mid point
                        %(for x, y, z separately) per 8 sub-cells at this level 
                        %(there are 8 sub-cells (not necessarily occupied) per 
                        %occupied cell at the previous level). The coordinates of 
                        %occupied cells at the previous level are stored in 
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
                    end                  
                    %Compute the 8 corner coordinates for each of the 8 sub-cells
                    %resulting from subdividing each of the occupied cells whose
                    %coordinates are stored in subcell_coords_all (if lvl_d > lvl +
                    %1) or in current_corner_coords (if lvl_d = lvl + 1)
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
                    end
                    %For each corner coordinate in subcell_coords, compute the
                    %control point associated with this corner, by interpolating 
                    %(tri-linear interpolation) between the original parent control
                    %points, current_ctrlpts. Then, for each set of 8 corner coordinates
                    %(i.e., for each sub-cell), check if all 8 of this cell's
                    %interpolated control points have the same sign: if they don't, store this
                    %sub-cell's coordinates in subcell_coords_all for the current
                    %level, "lvl_d"; if they do have the same sign, do nothing
                    %further, as this sub-cell will not be a candidate for further
                    %subdivision at the next octree level.
                    disp(['  Computing and checking the signs of interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
                    zcc_coords_cntr = 1;    %Counter for coordinates of cells at level "lvl_d", which contain zero crossings
                    for sub_cell_8set = 1:8:(size(subcell_coords, 1) - 7)
                        %Get all 8 corner coordinates for the current sub-cell
                        sc_coords = subcell_coords((sub_cell_8set:(sub_cell_8set + 7)), :);
                        %Normalize sc_coords to be in the range [0, 1], because the
                        %trilinear interpolation formula (used below) expects the
                        %(x, y, z) values to be in this range
                        sc_coords_orig = sc_coords;
                        %sc_coords = (sc_coords - sc_coords(1, :))/(2^(b + 1 - lvl_d));   %Subtract the origin of the current sub-cell, and divide by the cell width
                        sc_coords = (sc_coords - current_corner_coords(1, :))/(2^(b + 1 - lvl));
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
        %                 for c = 1:8
        %                     subcell_ctrlpts(c) = current_ctrlpts(1)*(1 - sc_coords(c, 1))*(1 - sc_coords(c, 2))*(1 - sc_coords(c, 3)) ...
        %                         + current_ctrlpts(2)*sc_coords(c, 1)*(1 - sc_coords(c, 2))*(1 - sc_coords(c, 3)) ...
        %                         + current_ctrlpts(3)*sc_coords(c, 1)*sc_coords(c, 2)*(1 - sc_coords(c, 3)) ...
        %                         + current_ctrlpts(4)*(1 - sc_coords(c, 1))*sc_coords(c, 2)*(1 - sc_coords(c, 3)) ...
        %                         + current_ctrlpts(5)*(1 - sc_coords(c, 1))*(1 - sc_coords(c, 2))*sc_coords(c, 3) ...
        %                         + current_ctrlpts(6)*sc_coords(c, 1)*(1 - sc_coords(c, 2))*sc_coords(c, 3) ...
        %                         + current_ctrlpts(7)*sc_coords(c, 1)*sc_coords(c, 2)*sc_coords(c, 3) ...
        %                         + current_ctrlpts(8)*(1 - sc_coords(c, 1))*sc_coords(c, 2)*sc_coords(c, 3);
        %                 end
                        %Check interpolated control point signs
                        if (sum(sign(subcell_ctrlpts)) == length(subcell_ctrlpts))||(sum(sign(subcell_ctrlpts)) == -length(subcell_ctrlpts))
                            %If all control points have the same sign
                            continue;
                        else
                            %If control points do not all have the same sign
                            subcell_coords_all{lvl_d, occ_cell}((zcc_coords_cntr:(zcc_coords_cntr + 7)), 1:3) = sc_coords_orig;
                            zcc_coords_cntr = zcc_coords_cntr + 8;
                        end
                    end      
                end %End lvl_d
            end %End check if lvl == (b + 1)  
        end %End check if control points for current occ_cell have same signs
        disp('------------------------------------------------------------');
    end %End occ_cell

    disp(['Collecting all reconstructed voxels for each occupied cell at level ' num2str(lvl) ' and finding their centre coordinates ...']);
    %Collect all of the sub-cell coordinates stored at the leaf levels of
    %subcell_coords_all, for each occ_cell at level "lvl": these represent
    %our reconstructed voxel (x, y, z) corner coordinates obtained from the 
    %reconstructed control points at level "lvl".
%     v_cntr = 1;
%     for sca = 1:size(subcell_coords_all, 2)
%         if ~isempty(subcell_coords_all{end, sca})
%             reconstructed_vox_pos_corners(v_cntr:(v_cntr + size(subcell_coords_all{end, sca}, 1) - 1), 1:3) = subcell_coords_all{end, sca};
%             v_cntr = v_cntr + size(subcell_coords_all{end, sca}, 1);
%         end
%     end
    reconstructed_vox_pos_corners = cat(1, subcell_coords_all{end, :});
    %Our input point cloud's voxel coordinates were considered to be the 
    %centres of the 1x1x1 voxels, so find the midpoint of each voxel cell
    %in reconstructed_vox_pos_corners: these midpoints will be our 
    %reconstructed voxel positions
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
    disp('------------------------------------------------------------');
    %Plot the reconstructed voxels
    figure;
    scatter3(reconstructed_vox_pos(:, 1), reconstructed_vox_pos(:, 2), reconstructed_vox_pos(:, 3), 5, 'filled');
    axis equal; axis off;
    title({'Voxel Reconstruction Using Only Reconstructed Control Points', ['at Octree Level ' num2str(lvl) ' and Interpolating at Higher Levels']});
    %Save the above reconstruction as a MATLAB figure and as a PDF image in
    %our network directory (NB: The '-bestfit' option maximizes the size of 
    %the figure to fill the page, but preserves the aspect ratio of the 
    %figure. The figure might not fill the entire page. This option leaves 
    %a minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_cplvl' num2str(lvl)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\vox_recon_cplvl' num2str(lvl)], '-dpdf');
    disp('Saving reconstructed voxels figure ...');
    disp('------------------------------------------------------------');
    disp(' ');
end
disp('------------------- DECODER FINISHED -----------------------');
disp(' ');













