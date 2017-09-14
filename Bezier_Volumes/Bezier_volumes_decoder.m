%Requires directory "phil": add this directory plus its sub-directories to
%the current MATLAB path.

%start_lvl = 1;
%max_lvl = 8;
%b = 10;
%rec_ctrlpts_start_lvl = reconstructed_control_points{start_lvl};
%OccupancyCode = myOT.OccupancyCode; %Should be from start_lvl only
%vis_levels_ot = 4; %No. of octree levels for which we want to visualize the octree cell subdivision

function [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, reconstructed_control_points, prune_flag, varargin)

disp(' ');
disp('============================================================');
disp('                   DECODER RUNNING ...');
disp('============================================================');
disp(' ');

OccupancyCode = occupancy_codes_forDec;
%post_pruning_array = post_pruning_array_forDec;

%Check if octree and wavelet coefficient tree pruning was used at the
%encoder
if numel(varargin) >= 1
    post_pruning_array = varargin{1};
end

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

if prune_flag == 1
    %Find the first non-empty location in post_pruning_array: this 
    %indicates the first octree level at which leaf cells are found
    %NOTE: Currently assuming that there are no octree levels AFTER
    %pp_fist_nonempty that do not contain any leaf cells.
    pp_first_nonempty = find(~cellfun(@isempty, post_pruning_array), 1);
end

%For each octree level ...
%for lvl = 1:b
for lvl = 1:(max_lvl - 1)
    disp(['Processing octree level ' num2str(lvl) ':']);
    %Counter for all children of all occupied nodes at this level
    total_child_cntr = 1;
    if ((prune_flag == 1)&&(lvl < pp_first_nonempty))||(prune_flag == 0)
        %For each occupied cell at this level ...
        for occ_cell = 1:numel(OccupancyCode{lvl})
            %Convert the OccupancyCode decimal value for this cell's 
            %children, into its binary representation
            bin_vec = dec2bin(OccupancyCode{lvl}(occ_cell), 8);
            %The number of "1"s in bin_vec indicates the number of occupied 
            %children that this octree cell has
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
                    %contains non-zero values: these are the positions to
                    %which we will add offsets to obtain the children's 
                    %spatial indices, when necessary
                    parent_SI_nonzero = uint16(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                    %If parent_SI_nonzero is empty (i.e., the current 
                    %occupied cell's spatial index is (0, 0, 0)), then the 
                    %spatial indices of this cell's children can be 
                    %obtained just by indexing into SI_dict
                    if isempty(parent_SI_nonzero)
                        SpatialIndex{lvl + 1}(1:ChildCount{lvl}(occ_cell), :) = uint16(SI_dict(ones_inds, :));
                    else
                        %Add an offset to corresponding locations in 
                        %SI_dict, in the columns determined by 
                        %parent_SI_nonzero
                        SI_dict_locations = SI_dict(ones_inds, :);
                        SI_dict_locations(:, parent_SI_nonzero) = uint16(SI_dict_locations(:, parent_SI_nonzero)) + uint16(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                        SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(occ_cell) - 1), :) = uint16(SI_dict_locations);
                    end
                end
            end %End occ_child  
            total_child_cntr = total_child_cntr + ChildCount{lvl}(occ_cell);
        end %End occ_cell
    elseif (prune_flag == 1)&&(lvl >= pp_first_nonempty)
        %Counter to keep track of where in the OccupancyCode array we are 
        %up to at the current octree level: the pruned OccupancyCode array 
        %only stores occupancy codes for internal (non-leaf) cells, so this
        %counter will only be incremented each time we come across an
        %internal cell
        oc_code_cntr = 0;
        %For each occupied cell at this level ...
        for occ_cell = 1:numel(post_pruning_array{lvl})
            %If this cell is a leaf
            if post_pruning_array{lvl}(occ_cell) == 1
                %It has no children (they were pruned off at the encoder),
                %so do nothing further
                continue;  
            %If this cell is not a leaf (i.e., it is internal)
            elseif post_pruning_array{lvl}(occ_cell) == 0
                %We advance 1 step in OccupancyCode array each time we find
                %an internal octree node
                oc_code_cntr = oc_code_cntr + 1;
                %Convert the OccupancyCode decimal value for this cell's 
                %children, into its binary representation
                bin_vec = dec2bin(OccupancyCode{lvl}(oc_code_cntr), 8);
                %The number of "1"s in bin_vec indicates the number of 
                %occupied children that this octree cell has
                ChildCount{lvl}(oc_code_cntr) = numel(strfind(bin_vec, '1'));  
                %Find the locations in bin_vec where the vector contains 
                %'1's
                ones_inds = strfind(bin_vec, '1');
                %Compute SpatialIndex for each occupied child    
                for occ_child = 1:length(ChildCount{lvl}(oc_code_cntr))
                    %If we are currently at the root cell, we can obtain 
                    %its children's spatial indices just by indexing into
                    %SI_dict
                    if lvl == 1
                        SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint16(SI_dict(ones_inds, :));
                    else
                        %Find where the curent occupied cell's spatial 
                        %index contains non-zero values: these are the 
                        %positions to which we will add offsets to obtain 
                        %the children's spatial indices, when necessary
                        parent_SI_nonzero = uint16(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                        %If parent_SI_nonzero is empty (i.e., the current 
                        %occupied cell's spatial index is (0, 0, 0)), then
                        %the spatial indices of this cell's children can be 
                        %obtained just by indexing into SI_dict
                        if isempty(parent_SI_nonzero)
                            SpatialIndex{lvl + 1}(1:ChildCount{lvl}(oc_code_cntr), :) = uint16(SI_dict(ones_inds, :));
                        else
                            %Add an offset to corresponding locations in 
                            %SI_dict, in the columns determined by 
                            %parent_SI_nonzero
                            SI_dict_locations = SI_dict(ones_inds, :);
                            SI_dict_locations(:, parent_SI_nonzero) = uint16(SI_dict_locations(:, parent_SI_nonzero)) + uint16(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                            SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(oc_code_cntr) - 1), :) = uint16(SI_dict_locations);
                        end
                    end
                end %End occ_child  
                total_child_cntr = total_child_cntr + ChildCount{lvl}(oc_code_cntr);
            end %End check if post_pruning_array{lvl}(occ_cell) == 1
        end %End occ_cell
    end %End check if lvl < pp_first_nonempty 
    %Make ChildCount{lvl} a column vector rather than a row vector
    ChildCount{lvl} = ChildCount{lvl}';
    disp('Finished computing ChildCount for each occupied cell at this level, and SpatialIndex for each occupied child of each occupied cell');
    %Compute the first child pointer for each occupied cell at the current
    %octree level
    lastChildPtr = cumsum(int32(ChildCount{lvl}));
    FirstChildPtr{lvl} = uint32(lastChildPtr - int32(ChildCount{lvl}) + 1);
    disp('Finished computing lastChildPtr and FirstChildPtr for each occupied cell at this level');
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

reconstruction_decoder = reconstruct_control_points_decoder(rec_ctrlpts_forDec, wavelet_coeffs_forDec, OccupancyCode, FirstChildPtr, ChildCount, corner_coords_decoder, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, prune_flag, reconstructed_control_points);

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

%-------- Point Cloud Voxel Reconstruction at Fixed Octree Levels --------%

%This can be done if no octree pruning was done at the encoder, because all
%the leaf voxels will be at the same octree level (b + 1), and so we can
%choose to use the reconstructed control points at certain octree level(s) 
%and interpolate between them for the higher octree levels, to figure 
%out through which octree cells the surface passes, right down to the voxel
%level.

if prune_flag == 0
    disp(' ');
    disp('-------- Voxel Reconstruction at a Fixed OT Level ----------');
    disp(' ');

    %[reconstructed_vox_pos, reconstructed_vox_pos_corners, subcell_coords_all] = voxel_reconstruction_nopruning(SpatialIndex, corner_coords_decoder, reconstruction_decoder, ctrl_pts_pointers, b, max_lvl, ptcloud_file, ptcloud_name);
    [reconstructed_vox_pos, ~, ~] = voxel_reconstruction_nopruning(SpatialIndex, corner_coords_decoder, reconstruction_decoder, ctrl_pts_pointers, b, max_lvl, ptcloud_file, ptcloud_name);
end

%--------------- Voxel Reconstruction using Pruned Octree ----------------%

%In a pruned octree, the leaf cells can be at different octree levels, so
%the below code will generate only one final point cloud reconstruction,
%using all of the control points that were sent to the decoder and
%interpolating only where the leaf cells are not already at the voxel level 
%(i.e., they are at level < b + 1).

if prune_flag == 1
    disp(' ');
    disp('-------- Voxel Reconstruction using Pruned Octree ----------');
    disp(' ');
    
    %[reconstructed_vox_pos, reconstructed_vox_pos_corners, subcell_coords_all] = voxel_reconstruction_pruning(pp_first_nonempty, SpatialIndex, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, ptcloud_file, ptcloud_name);
    [reconstructed_vox_pos, ~] = voxel_reconstruction_pruning(pp_first_nonempty, SpatialIndex, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, ptcloud_file, ptcloud_name);
end

disp(' ');
disp('------------------- DECODER FINISHED -----------------------');
disp(' ');













