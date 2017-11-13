%Requires directory "phil": add this directory plus its sub-directories to
%the current MATLAB path.

%start_lvl = 1;
%max_lvl = 8;
%b = 10;
%rec_ctrlpts_start_lvl = reconstructed_control_points{start_lvl};
%OccupancyCode = myOT.OccupancyCode; %Should be from start_lvl only
%vis_levels_ot = 4; %No. of octree levels for which we want to visualize the octree cell subdivision

function [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(debug_flag, occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, reconstructed_control_points, prune_flag, varargin)

disp(' ');
disp('============================================================');
disp('                   DECODER RUNNING ...');
disp('============================================================');
disp(' ');

start_dec_time = tic;

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

start_OT_recon_time = tic;

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
%SpatialIndex{1} = [uint16(0), uint16(0), uint16(0)];
SpatialIndex{1} = [uint32(0), uint32(0), uint32(0)];

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
    %indicates the first octree level at which leaf cells are found.
    %NOTE: Currently assuming that there are no octree levels AFTER
    %pp_fist_nonempty that do not contain any leaf cells.
    pp_first_nonempty = find(~cellfun(@isempty, post_pruning_array), 1);
end

%For each octree level ...
%for lvl = 1:b
for lvl = 1:(max_lvl - 1)
    if debug_flag == 1
        disp(['Processing octree level ' num2str(lvl) ':']);
    end
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
                    %SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint16(SI_dict(ones_inds, :));
                    SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint32(SI_dict(ones_inds, :));
                else
                    %Find where the curent occupied cell's spatial index 
                    %contains non-zero values: these are the positions to
                    %which we will add offsets to obtain the children's 
                    %spatial indices, when necessary
                    %parent_SI_nonzero = uint16(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                    parent_SI_nonzero = uint32(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                    %If parent_SI_nonzero is empty (i.e., the current 
                    %occupied cell's spatial index is (0, 0, 0)), then the 
                    %spatial indices of this cell's children can be 
                    %obtained just by indexing into SI_dict
                    if isempty(parent_SI_nonzero)
                        %SpatialIndex{lvl + 1}(1:ChildCount{lvl}(occ_cell), :) = uint16(SI_dict(ones_inds, :));
                        SpatialIndex{lvl + 1}(1:ChildCount{lvl}(occ_cell), :) = uint32(SI_dict(ones_inds, :));
                    else
                        %Add an offset to corresponding locations in 
                        %SI_dict, in the columns determined by 
                        %parent_SI_nonzero
                        SI_dict_locations = SI_dict(ones_inds, :);
                        %SI_dict_locations(:, parent_SI_nonzero) = uint16(SI_dict_locations(:, parent_SI_nonzero)) + uint16(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                        SI_dict_locations(:, parent_SI_nonzero) = uint32(SI_dict_locations(:, parent_SI_nonzero)) + uint32(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                        %SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(occ_cell) - 1), :) = uint16(SI_dict_locations);
                        SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(occ_cell) - 1), :) = uint32(SI_dict_locations);
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
            end
            %If this cell is not a leaf (i.e., it is internal), we advance
            %1 step in the OccupancyCode array
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
                    %SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint16(SI_dict(ones_inds, :));
                    SpatialIndex{lvl + 1}(1:ChildCount{1}, :) = uint32(SI_dict(ones_inds, :));
                else
                    %Find where the curent occupied cell's spatial 
                    %index contains non-zero values: these are the 
                    %positions to which we will add offsets to obtain 
                    %the children's spatial indices, when necessary
                    %parent_SI_nonzero = uint16(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                    parent_SI_nonzero = uint32(find(SpatialIndex{lvl}(occ_cell, :) > 0));
                    %If parent_SI_nonzero is empty (i.e., the current 
                    %occupied cell's spatial index is (0, 0, 0)), then
                    %the spatial indices of this cell's children can be 
                    %obtained just by indexing into SI_dict
                    if isempty(parent_SI_nonzero)
                        %SpatialIndex{lvl + 1}(1:ChildCount{lvl}(oc_code_cntr), :) = uint16(SI_dict(ones_inds, :));
                        SpatialIndex{lvl + 1}(1:ChildCount{lvl}(oc_code_cntr), :) = uint32(SI_dict(ones_inds, :));
                    else
                        %Add an offset to corresponding locations in 
                        %SI_dict, in the columns determined by 
                        %parent_SI_nonzero
                        SI_dict_locations = SI_dict(ones_inds, :);
                        %SI_dict_locations(:, parent_SI_nonzero) = uint16(SI_dict_locations(:, parent_SI_nonzero)) + uint16(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                        SI_dict_locations(:, parent_SI_nonzero) = uint32(SI_dict_locations(:, parent_SI_nonzero)) + uint32(2*SpatialIndex{lvl}(occ_cell, parent_SI_nonzero));
                        %SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(oc_code_cntr) - 1), :) = uint16(SI_dict_locations);
                        SpatialIndex{lvl + 1}(total_child_cntr:(total_child_cntr + ChildCount{lvl}(oc_code_cntr) - 1), :) = uint32(SI_dict_locations);
                    end
                end
            end %End occ_child  
            total_child_cntr = total_child_cntr + ChildCount{lvl}(oc_code_cntr);
        end %End occ_cell
    end %End check if lvl < pp_first_nonempty 
    %Make ChildCount{lvl} a column vector rather than a row vector
    ChildCount{lvl} = ChildCount{lvl}';
    if debug_flag == 1
        disp('Finished computing ChildCount for each occupied cell at this level, and SpatialIndex for each occupied child of each occupied cell');
    end
    %Compute the first child pointer for each occupied cell at the current
    %octree level
    lastChildPtr = cumsum(int32(ChildCount{lvl}));
    FirstChildPtr{lvl} = uint32(lastChildPtr - int32(ChildCount{lvl}) + 1);
    if debug_flag == 1
        disp('Finished computing lastChildPtr and FirstChildPtr for each occupied cell at this level');
        disp('------------------------------------------------------------');
    end
end %End lvl
OT_recon_time = toc(start_OT_recon_time);
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
%Find the first octree level at which SpatialIndex is empty (assume that
%SpatialIndex is also empty at all levels higher than first_SI_empty)
first_SI_empty = find(cellfun(@isempty, SpatialIndex), 1);

%Figure out the corner coordinates for each corner of each occupied cell at
%each octree level, starting from start_lvl
start_cornercoords_time = tic;
%for lvl = start_lvl:(b + 1)
%for lvl = start_lvl:max_lvl
if isempty(first_SI_empty)
    end_lvl = max_lvl;
else
    end_lvl = first_SI_empty - 1;
end
for lvl = start_lvl:end_lvl
    if debug_flag == 1
        disp(['Processing octree level ' num2str(lvl) ':']); 
        disp('Computing corner coordinates for each occupied cell ...');
        disp('------------------------------------------------------------');
    end
    
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
cornercoords_time = toc(start_cornercoords_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all corner coordinates: ' num2str(cornercoords_time) ' seconds']);
disp('************************************************************');

%---------------- Signal (Control Point) Reconstruction ------------------%

disp(' ');
disp('-------------- Control Point Reconstruction ----------------');
disp(' ');

if prune_flag == 0
    reconstruction_decoder = reconstruct_control_points_decoder(debug_flag, rec_ctrlpts_forDec, wavelet_coeffs_forDec, SpatialIndex, FirstChildPtr, ChildCount, corner_coords_decoder, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, prune_flag, reconstructed_control_points);
elseif prune_flag == 1
    reconstruction_decoder = reconstruct_control_points_decoder(debug_flag, rec_ctrlpts_forDec, wavelet_coeffs_forDec, SpatialIndex, FirstChildPtr, ChildCount, corner_coords_decoder, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, prune_flag, reconstructed_control_points, post_pruning_array, pp_first_nonempty);
end

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
    [reconstructed_vox_pos, ~] = voxel_reconstruction_pruning(debug_flag, pp_first_nonempty, corner_coords_decoder, post_pruning_array, reconstruction_decoder, ctrl_pts_pointers, b, q_stepsize);
end

total_decoder_time = toc(start_dec_time);

disp(' ');
disp('----------- Displaying Reconstruction Results --------------');
disp(' ');

%Plot the reconstructed voxels
figure;
%NOTE: Below, we are only reading in the original PLY file in order to get 
%the corresponding colours assigned to each reconstructed voxel (this will
%only work when the same number of voxels is reconstructed as in the 
%original point cloud, and these voxels are in the same order as the 
%original voxels, which will be done below)
[~, ptcloud, ~] = plyRead(ptcloud_file);
if size(reconstructed_vox_pos, 1) == size(ptcloud, 1)
    %Order the reconstructed voxels according to their Morton codes, so
    %that they are in the same order as the input point cloud at the
    %encoder
    if debug_flag == 1
        disp('Reordering reconstructed voxels according to Morton codes ...');
    end
    %Get Morton codes for the reconstructed voxel x, y, z coordinates
    mortonCodes = xyzToMorton(reconstructed_vox_pos, b);   %"b" bits for each Morton code
    if debug_flag == 1
        disp('Morton codes computed');
    end
    %Sort the Morton codes obtained above, in ascending order
    [~, I_vox] = sort(mortonCodes);
    if debug_flag == 1
        disp('Morton codes sorted');
    end
    %Sort the voxel x, y, z locations in the same order as the sorted 
    %Morton codes
    reconstructed_vox_pos = reconstructed_vox_pos(I_vox, 1:3);
    if debug_flag == 1
        disp('Reconstructed voxels sorted');
        disp('------------------------------------------------------------');
    end
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
    title({'Original Voxels that were NOT Reconstructed', ['at Decoder (' num2str(size(test_vox_diffs, 1)) '/' num2str(size(ptcloud, 1)) ' = ' num2str((size(test_vox_diffs, 1)/size(ptcloud, 1))*100) '%)']});
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
disp(['Total number of reconstructed voxels NOT present in the original point cloud: ' num2str(size(test_vox_diffs2, 1))]);
if ~isempty(test_vox_diffs2)
    figure;
    scatter3(test_vox_diffs2(:, 1), test_vox_diffs2(:, 2), test_vox_diffs2(:, 3), 5, 'filled', 'MarkerFaceColor', 'r');
    axis equal; axis off;
    title(['Reconstructed Voxels that were NOT in Input Point Cloud: ' num2str(size(test_vox_diffs2, 1))]);
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_post_pruning']);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\incorrect_voxels_post_pruning'], '-dpdf');
    disp('Saving surplus/incorrect voxels figure ...');
    disp('------------------------------------------------------------');
end   

disp(' ');
disp('------------------- DECODER FINISHED -----------------------');
disp(' ');

disp(['TOTAL DECODER TIME (excluding displaying of results): ' num2str(total_decoder_time) ' seconds']);
disp(' ');













