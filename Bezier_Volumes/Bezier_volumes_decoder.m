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
    %indicates the first octree level at which octree pruning occurred at 
    %the encoder
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

%------ Point Cloud Voxel Reconstruction at Different Octree Levels ------%

%This can be done if no octree pruning was done at the encoder, because all
%the leaf voxels will be at the same octree level (b + 1).

if prune_flag == 0
    disp(' ');
    disp('-------- Voxel Reconstruction at a Fixed OT Level ----------');
    disp(' ');

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
end %End check if prune_flag == 0

%--------------- Voxel Reconstruction using Pruned Octree ----------------%

%In a pruned octree, the leaf cells can be at different octree levels, so
%the below code will generate only one final point cloud reconstruction,
%using all of the control points that were sent to the decoder and
%interpolating only where the leaf cells are not already at the voxel level 
%(i.e., they are at level b < b + 1).

if prune_flag == 1
    disp(' ');
    disp('-------- Voxel Reconstruction using Pruned Octree ----------');
    disp(' ');
    
    %Cell array to store the coordinates of octree cells at all levels,  
    %which are candidates for further subdivision
    subcell_coords_all = cell((b + 1), size(SpatialIndex{lvl}, 1)); 
    %Matrix to store the reconstructed voxel corners
    reconstructed_vox_pos_corners = [];
%     %Counter to keep track of the number of reconstructed voxel corners
%     rec_vox_cnr_cntr = 1;
    
    %For each octree level, starting from the level where the first leaf
    %cells might be found due to pruning ...
    for lvl = pp_first_nonempty:size(SpatialIndex, 1)
        disp(' ');
        disp(['Reconstructing voxels for level ' num2str(lvl) ' leaf cells ...']);
        disp(' ');
        %If we are at the voxel level ...
        if lvl == (b + 1)
            %Just get all of the occupied voxels' corner coordinates and
            %store them  
            disp(['Retrieving reconstructed corner coordinates for occupied voxels at level ' num2str(lvl)]);
            reconstructed_vox_pos_corners = corner_coords_decoder{lvl};
            continue;
        end

        %Counter for the number of leaf octree cells found at the current
        %level (for debugging purposes only)
        nbr_leaves = 0;
        %For each occupied cell at the current level ...
        for occ_cell = 1:size(SpatialIndex{lvl}, 1)
%             %If we are at the final leaf level (i.e., voxel level)
%             if lvl == (b + 1)
%                 nbr_leaves = nbr_leaves + 1;
%                 %Just get the current occupied voxel's corner coordinates 
%                 %and store them  
%                 %disp(['Retrieving corner coordinates for occupied voxel ' num2str(occ_cell)]);
%                 current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
%                 reconstructed_vox_pos_corners((rec_vox_cnr_cntr:(rec_vox_cnr_cntr + 7)), 1:3) = current_corner_coords;
%                 rec_vox_cnr_cntr = rec_vox_cnr_cntr + 8;
%                 continue;
%             %If we are not yet at the final leaf level (voxel level)
%             else
                %If the current cell is not a leaf
                if post_pruning_array{lvl}(occ_cell) == 0
                    %Keep going until we find a leaf
                    continue;
                %If the current cell is a leaf
                elseif post_pruning_array{lvl}(occ_cell) == 1
                    nbr_leaves = nbr_leaves + 1;
                    disp(['Reconstructing voxels for leaf cell (occ_cell) ' num2str(occ_cell) ':']);
                    %Get the corner coordinates for each of the 8 corners 
                    %of this cell
                    current_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :);
                    %Get the control points for all 8 corners of the 
                    %current cell ...
                    current_ctrlpts = reconstruction_decoder{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
                    %Subdivide the current occ_cell until we get to a cell 
                    %of size 1x1x1 (i.e., a voxel), and interpolate between 
                    %the original cell's (occ_cell's) control points at 
                    %each successive subdivision to figure out which voxels 
                    %are occupied at the end ... 
                    %For each descendant of the current occ_cell, right 
                    %down to the final leaf level (voxel level) ...
                    for lvl_d = (lvl + 1):(b + 1)   
                        disp(['--Processing sub-cells at level ' num2str(lvl_d)]);
                        %Find the unique minimum, unique maximum, and 
                        %midpoint coordinates for the x, y, and z 
                        %dimensions of the sub-cells at this level
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
                            %There will be one unq_min, one unq_max, and 
                            %one mid point (for x, y, z separately) per 8 
                            %sub-cells at this level (there are 8 sub-cells 
                            %(not all necessarily occupied) per occupied 
                            %cell at the previous level). The coordinates 
                            %of cells at the previous level, which are 
                            %candidates for further subdivision, are stored 
                            %in subcell_coords_all(lvl_d - 1).
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
                        %sub-cells resulting from subdividing each of the 
                        %cells whose coordinates are stored in 
                        %subcell_coords_all (if lvl_d > lvl + 1) or in 
                        %current_corner_coords (if lvl_d = lvl + 1)
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
                        %For each corner coordinate in subcell_coords, 
                        %compute the control point associated with this 
                        %corner by interpolating (tri-linearly) between 
                        %the original parent control points, 
                        %current_ctrlpts
                        disp(['  Computing interpolated control points for each sub-cell at this level (8 control points per sub-cell, ' num2str(size(subcell_coords, 1)/8) ' sub-cells in total)']);
                        zcc_coords_cntr = 1;    %Counter for coordinates of cells at lvl_d, which will be candidates for further subdivision
                        for sub_cell_8set = 1:8:(size(subcell_coords, 1) - 7)
                            %Get all 8 corner coordinates for the current 
                            %sub-cell
                            sc_coords = subcell_coords((sub_cell_8set:(sub_cell_8set + 7)), :);
                            %Normalize sc_coords to be in the range [0, 1],  
                            %because the trilinear interpolation formula 
                            %(used below) expects the (x, y, z) values to 
                            %be in this range
                            sc_coords_orig = sc_coords;
                            %sc_coords = (sc_coords - sc_coords(1, :))/(2^(b + 1 - lvl_d));   %Subtract the origin of the current sub-cell, and divide by the cell width
                            sc_coords = (sc_coords - current_corner_coords(1, :))/(2^(b + 1 - lvl));
                            %Initialize an array to store the interpolated 
                            %control points for all 8 corners of the 
                            %current sub-cell
                            subcell_ctrlpts = zeros(8, 1);
                            %Compute all 8 control points for the corners
                            %of the current sub-cell
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

                            %If we are one level before the leaf level, or
                            %at the leaf level
                            if lvl_d >= b
                                %Check interpolated control point signs
                                if (sum(sign(subcell_ctrlpts)) == length(subcell_ctrlpts))||(sum(sign(subcell_ctrlpts)) == -length(subcell_ctrlpts))
                                    %If all the control points have the 
                                    %same sign, consider the current 
                                    %sub-cell unoccupied: it will not be 
                                    %subdivided further
                                    continue;
                                else
                                    %If the control points do not all have
                                    %the same sign, consider the current 
                                    %sub-cell occupied: it will be 
                                    %subdivided further unless we are at 
                                    %the leaf level (b + 1)
                                    subcell_coords_all{lvl_d, occ_cell}((zcc_coords_cntr:(zcc_coords_cntr + 7)), 1:3) = sc_coords_orig;
                                    zcc_coords_cntr = zcc_coords_cntr + 8;
                                end
                            %If we are not yet at the leaf level or at one
                            %level before the leaf level
                            else
                                %The current sub-cell will be subdivided 
                                %further, regardless of its control point 
                                %signs
                                subcell_coords_all{lvl_d, occ_cell}((zcc_coords_cntr:(zcc_coords_cntr + 7)), 1:3) = sc_coords_orig;
                                zcc_coords_cntr = zcc_coords_cntr + 8;
                            end
                        end %End sub_cell_8set      
                    end %End lvl_d
                end %End check if post_pruning_array{lvl}(occ_cell) == 1
%             end %End check if lvl == (b + 1)
        end %End occ_cell
        disp('------------------------------------------------------------');
        disp(['Total number of leaf cells at this level: ' num2str(nbr_leaves)]);
        disp('------------------------------------------------------------');
    end %End lvl
    
    %Collect all of the sub-cell coordinates stored at level b of 
    %subcell_coords_all, for each occ_cell: these represent our 
    %reconstructed voxel (x, y, z) corner coordinates obtained from the 
    %reconstructed control points
    disp('Collecting all reconstructed voxels for each occupied cell and finding their centre coordinates ...');
    reconstructed_vox_pos_corners2 = cat(1, subcell_coords_all{end, :});
    %Concatenate the reconstructed voxel corners already stored in
    %reconstructed_vox_pos with the other corners collected above
    reconstructed_vox_pos_corners = cat(1, reconstructed_vox_pos_corners, reconstructed_vox_pos_corners2);
    
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
        %separately): these mean values represent the centre (x, y, z)
        %location of the current voxel
        reconstructed_vox_pos(vox_cntr, :) = mean(vc_coords, 1);
        vox_cntr = vox_cntr + 1;
    end

%     %Order the reconstructed voxels according to their Morton codes, so
%     %that they are in the same order as the input point cloud at the
%     %encoder
%     disp('Reordering reconstructed voxels according to Morton codes ...');
%     %Get Morton codes for the reconstructed voxel x, y, z coordinates
%     mortonCodes = xyzToMorton(reconstructed_vox_pos, b);   %"b" bits for each Morton code
%     disp('Morton codes computed');
%     %Sort the Morton codes obtained above, in ascending order
%     [~, I_vox] = sort(mortonCodes);
%     disp('Morton codes sorted');
%     %Sort the voxel x, y, z locations in the same order as the sorted 
%     %Morton codes
%     reconstructed_vox_pos = reconstructed_vox_pos(I_vox, 1:3);
%     disp('Reconstructed voxels sorted');
%     disp('------------------------------------------------------------');
    
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

end %End check if prune_flag == 0

disp(' ');
disp('------------------- DECODER FINISHED -----------------------');
disp(' ');













