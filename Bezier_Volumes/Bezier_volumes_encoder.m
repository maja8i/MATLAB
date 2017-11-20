%DOES WAVELET ANALYSIS FROM THE ROOT LEVEL TO (LEAVES-1), QUANTIZING ALONG
%THE WAY.

%Requires directory "phil": add this directory plus its sub-directories to
%the current MATLAB path.

% ptcloud_file = '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10_WithNormals\boxer_voxelized10.ply';
% b = 10;
% start_lvl = 1;
% max_lvl = 8;
% vis_levels_ot = 5; %No. of octree levels for which we want to visualize the octree cell subdivision
% vis_levels_ctrlpts = 2; %No. of octree levels for which we want to visualize the control point computation

function [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points, varargout] = Bezier_volumes_encoder(debug_flag, ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag, varargin)

disp(' ');
disp('============================================================');
disp('                   ENCODER RUNNING ...');
disp('============================================================');
disp(' ');

%Read in the input point cloud (assume PLY format)
[~, ptcloud, ~] = plyRead(ptcloud_file);

start_enc_time = tic;

if ~isempty(varargin)
    %Threshold for pruning wavelet coefficients
    zero_threshold = varargin{1};
end

%-------------------------- Octree Construction --------------------------%

disp('-------------------- Octree Construction -------------------');
disp(' ');

%Construct an octree of depth (b + 1), for the input point cloud
[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted] = construct_octree(debug_flag, ptcloud_file, b, ptcloud_name);
%[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted] = construct_octree(ptcloud_file, b);

disp(' ');
disp('-------------- Computing Corner Coordinates ----------------');
disp(' ');

%Initialize a cell array to store the corner coordinates of each occupied
%cell at each level of the octree
corner_coords = cell((b + 1), 1);
%Initialize a cell array to store only the unique corner coordinates at
%each octree level
unique_coords = cell((b + 1), 1);
%Initialize a cell array of "pointers": there will be 8 pointers per
%occupied octree cell (1 per corner) at each octree level, which will point
%to the Bezier control point in control_points (computed in 
%compute_control_points(), later), associated with that corner. So, for 
%each octree level "lvl", there will be myOT.NodeCount(lvl) pointer arrays
%containining 8 elements each.  
ctrl_pts_pointers = cell((b + 1), 1);

start_cornercoords_time = tic;
%For each octree level ...
%for lvl = 1:(b + 1) 
for lvl = 1:max_lvl 
    if debug_flag == 1
        disp(['Processing octree level ' num2str(lvl) ':']); 
        disp('Computing corner coordinates for each occupied cell ...');
        disp('------------------------------------------------------------');
    end
    
    %Find the (x, y, z) coordinates of the origin of each occupied octree
    %cell at the current level
    corners1 = double(myOT.SpatialIndex{lvl})*(2^(b + 1 - lvl)) - [0.5 0.5 0.5];
    %Replicate each row of corners1 7 times, so that we can directly add
    %this matrix to offsets_from_origin_rep (below)
    corners1_rep = corners1(repmat(1:size(corners1, 1), 7, 1), :);
    
    %Find the (x, y, z) coordinates of the other 7 corners of each occupied
    %octree cell at the current level
    offsets_from_origin = [[2^(b + 1 - lvl) 0 0]; [2^(b + 1 - lvl) 2^(b + 1 - lvl) 0]; [0 2^(b + 1 - lvl) 0]; [0 0 2^(b + 1 - lvl)]; [2^(b + 1 - lvl) 0 2^(b + 1 - lvl)]; [2^(b + 1 - lvl) 2^(b + 1 - lvl) 2^(b + 1 - lvl)]; [0 2^(b + 1 - lvl) 2^(b + 1 - lvl)]];
    offsets_from_origin_rep = repmat(offsets_from_origin, size(corners1, 1), 1);
    corners2_8 = corners1_rep + offsets_from_origin_rep;
    
    %Store all the corner coordinates for the current level in their
    %corresponding locations inside corner_coords (want all 8 corner
    %coordinates for one octree cell to be listed one after the other, in
    %an 8 x 3 matrix inside corner_coords{lvl})
    corner_coords{lvl} = zeros((myOT.NodeCount(lvl)*8), 3);
    corner_coords{lvl}(1:8:(myOT.NodeCount(lvl)*8 - 7), :) = corners1;   
    corners2_8_cntr = 1;
    for next_ind = 2:8:(myOT.NodeCount(lvl)*8 - 6)
        corner_coords{lvl}((next_ind:(next_ind + 6)), :) = corners2_8((corners2_8_cntr:(corners2_8_cntr + 6)), :);
        corners2_8_cntr = corners2_8_cntr + 7;
    end

    %Find all the unique corner coordinates at the current level of the
    %octree, and in ctrl_pts_pointers{lvl}, for each row that represents an
    %(x, y, z) coordinate from corner_coords{lvl}, store a pointer for this
    %coordinate into the unique_coords{lvl} array, to say what the
    %coordinate triplet in this row should be. These pointers will also act
    %as pointers to the control points stored in control_points{lvl}, since
    %control_points{lvl} will be the same size as unique_coords{lvl}, as
    %each unique corner coordinate triplet will have one corresponding 
    %Bezier control point (signed distance value) associated with it.
    [unique_coords{lvl}, ~, ctrl_pts_pointers{lvl}] = unique(corner_coords{lvl}, 'rows', 'stable');
end %end "lvl"
cornercoords_time = toc(start_cornercoords_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all corner coordinates: ' num2str(cornercoords_time) ' seconds']);
disp('************************************************************');

% %--------------- Visualization of Octree Cell Subdivision ----------------%
% 
% vis_levels_ot = max_lvl;
% 
% disp(' ');
% disp('Displaying octree subdivision ...'); 
% disp('------------------------------------------------------------');
% 
% %Read in and display the input point cloud (assume PLY format)
% [~, ptcloud, ~] = plyRead(ptcloud_file);
% figure;
% scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
% axis equal; axis off;
% hold on;
% 
% %Define a colour matrix for each octree level for which we wish to 
% %visualize the octree cell subdivision
% colours = hsv(vis_levels_ot);   %Each row contains 3 columns (R, G, B)
% 
% %Initialize a cell array to store, for each unique corner vertex at each
% %octree level, a list of indices of the corner vertices that share an edge
% %with that vertex. 
% shared_edge_inds = cell(b, 1);
% 
% %Initialize an array to store handles for each octree level plot, to use
% %for the legend 
% h_ot = zeros(1, vis_levels_ot);
% 
% %Create a cell array of labels for the plotted octree levels, to use in the
% %legend
% otplot_legend = cell(1, vis_levels_ot);
% for l = 1:vis_levels_ot
%     if l == 1
%         otplot_legend{l} = 'Level 1 (Root)';
%     else
%         otplot_legend{l} = ['Level ' num2str(l)];
%     end
% end
% 
% %for lvl = 1:vis_levels_ot
% for lvl = [4, 5]
%     %Plot all the unique corner points at the current octree level, over
%     %the top of the input 3D point cloud
%     scatter3(unique_coords{lvl}(:, 1), unique_coords{lvl}(:, 2), unique_coords{lvl}(:, 3), 10, colours(lvl, :), 'filled');
%     hold on;
%     For two of the unique vertices to share an edge, they must have at
%     least one coordinate (x, or y, or z) in common. Work out all the edges
%     for the unique vertices, and use them to display the horizonal and 
%     vertical octree grid lines that show the cell divisions at the current
%     octree level (i.e., connect the points that share each edge with a
%     straight line).
%     for v = 1:size(unique_coords{lvl}, 1)
%         %Initialize an array to store the indices of the shared vertices,
%         %for the current vertex
%         edge_pts_indices = [];
%         %Initialize a counter to keep track of how many vertex indices
%         %corresponding to shared edges have been found
%         edge_pts_indices_cntr = 1;
%         %Get the coordinates of the current vertex
%         current_vtx = unique_coords{lvl}(v, :);
%         %Find all the vertices (amongst the list of unique corner vertices
%         %at the current octree level) that share an edge with current_vtx
%         [x_same_rows, ~] = find(current_vtx(1) == unique_coords{lvl}(:, 1));
%         x_same_rows(find(x_same_rows == v)) = [];   %Remove the current vertex index
%         [y_same_rows, ~] = find(current_vtx(2) == unique_coords{lvl}(:, 2));
%         y_same_rows(find(y_same_rows == v)) = [];   %Remove the current vertex index
%         [z_same_rows, ~] = find(current_vtx(3) == unique_coords{lvl}(:, 3));
%         z_same_rows(find(z_same_rows == v)) = [];   %Remove the current vertex index
%         temp_cat = [x_same_rows; y_same_rows; z_same_rows];
%         for v_ind = 1:length(temp_cat)
%             %Only keep the vertex indices that appear twice in temp_cat
%             if length(find(temp_cat == temp_cat(v_ind))) == 2
%                 edge_pts_indices(edge_pts_indices_cntr) = temp_cat(v_ind);
%                 edge_pts_indices_cntr = edge_pts_indices_cntr + 1;
%             end
%         end
%         %Keep only the unique indices in edge_pts_indices
%         edge_pts_indices = (unique(edge_pts_indices))';
%         %Get the (x, y, z) coordinates corresponding to all the
%         %edge_pts_indices
%         edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%         %Connect the current_vtx to all the edge_pts_coords with straight
%         %lines
%         for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%             h_ot(lvl) = plot3([current_vtx(1) edge_pts_coords(nv, 1)], [current_vtx(2) edge_pts_coords(nv, 2)], [current_vtx(3) edge_pts_coords(nv, 3)], 'Color', colours(lvl, :));
%             hold on;
%         end
%         %Store the current list of edge indices inside shared_edge_inds,
%         %for future reference
%         shared_edge_inds{lvl, v} = edge_pts_indices;
%     end
% end
% hold off;
% legend(h_ot, otplot_legend, 'Location', 'best');
% title({'Octree Subdivision and Computing Corner Coordinates of Octree Cells', '(ENCODER)'}, 'Interpreter', 'none');

%-------------------- Computing Bezier Control Points --------------------%

disp(' ');
disp('-------------- Extracting Voxel Coordinates ----------------');
disp(' ');

%Extract the set of occupied voxel coordinates (x, y, z) at all levels of
%the octree myOT 
[~, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_normal_averages, occupied_voxel_centroids, occupied_voxel_centroid_averages] = extract_occupied_voxels(debug_flag, myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted);
%[~, occupied_voxel_coords, occupied_voxel_normals] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, normals_sorted);

disp(' ');
disp('------------ Computing Bezier Control Points ---------------');
disp(' ');

[control_points, nearest_voxels, normal_nearest_vox, min_euclid_dist, difference_vectors, dot_products, thresh] = compute_control_points(debug_flag, myOT, corner_coords, unique_coords, ctrl_pts_pointers, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_centroids, occupied_voxel_normal_averages, occupied_voxel_centroid_averages, b, max_lvl);
%[control_points, ~, ~, ~, ~, ~, thresh] = compute_control_points(myOT, corner_coords, unique_coords, ctrl_pts_pointers, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_centroids, occupied_voxel_normal_averages, occupied_voxel_centroid_averages, b, max_lvl);

if debug_flag == 1
    disp(' ');
    %For debugging purposes, check how many of the control points (out of 
    %the control points for the UNIQUE occupied cell corners only) at each 
    %octree level, if any, are 0
    for lvl = start_lvl:1:max_lvl
        zero_ctrlpt_cnt = length(find(control_points{lvl} == 0));
        disp(['No. of 0 control points at level ' num2str(lvl) ', before quantization: ' num2str(zero_ctrlpt_cnt) '/' num2str(length(control_points{lvl}))]);
    end 
    disp(' ');

    %For debugging purposes, print out how many occupied octree cells at 
    %each octree level, if any, end up having all 8 of their control points 
    %with the same sign (+/-/0) ...
    %[~, ptcloud, ~] = plyRead(ptcloud_file);
    %Accumulate the control points for ALL the occupied octree cell corners 
    %at each level (not just the control points for the unique corners) 
    %into one cell array
    all_ctrlpts = get_all_ctrlpts(control_points, ctrl_pts_pointers, start_lvl, max_lvl); 
    same_sign_voxels = [];
%     cells = [];
    for lvl = start_lvl:max_lvl
        same_sign_cntr = 0;
        cell_cntr = 0;
        zero_cp_cntr = 0;
        for i = 1:8:(length(all_ctrlpts{lvl}) - 7) 
            cell_cntr = cell_cntr + 1;
            %Get all 8 control points for the corners of the current cell
            current_ctrlpts = all_ctrlpts{lvl}(i:(i + 7));
            %Check if all control points of the current cell have the same
            %sign, including the case where all the control points may be 0
            if (abs(sum(sign(current_ctrlpts))) == 8)||(~any(sign(current_ctrlpts)))
                same_sign_cntr = same_sign_cntr + 1;
                %disp(['Cell ' num2str(cell_cntr) ' has all control points with the same sign: ']);
                if ~any(sign(current_ctrlpts))
                    zero_cp_cntr = zero_cp_cntr + 1;
                end
                if lvl == b + 1
                    %Store this voxel, for debugging purposes
                    same_sign_voxels((size(same_sign_voxels, 1) + 1), 1:3) = ptcloud(cell_cntr, 1:3);
                end
    %             %Display the control points for all corners of this cell
    %             disp(num2str(current_ctrlpts));
    %             disp(' ');
    %             %Get the corner coordinates for each corner of this cell
    %             cell_corners = corner_coords{lvl}((i:(i + lvl - 1)), 1:3);
    %             disp('Cell corner coordinates:');
    %             disp(num2str(cell_corners));
    %             disp(' ');
        %         %Get the coordinates (either midpoint or centroid, whichever one of
        %         %these happened to be used in the control point computation) of the
        %         %nearest voxel found for each corner of the current voxel
        %         curr_nearest_vox = nearest_voxels{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
        %         disp('Nearest voxel found for each corner of the current voxel: ');
        %         disp(num2str(curr_nearest_vox));
        %         disp(' ');
        %         %Get the normal vector for the current voxel, from the input
        %         %voxelized point cloud
        %         disp('Voxel normal:');
        %         disp(normals_sorted(vox, :));
        %         %Get the normal vector for each nearest voxel found above
        %         curr_normal_vec = normal_nearest_vox{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
        %         disp('Nearest voxel normal for each corner: ');
        %         disp(num2str(curr_normal_vec));
        %         disp(' ');
        %         %Get the difference vector for each of the 8 corners of this voxel
        %         %(difference vector between each corner and the chosen voxel)
        %         vox_diff_vecs = difference_vectors{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
        %         disp('Voxel corner difference vectors:');
        %         disp(num2str(vox_diff_vecs));
        %         disp(' ');
        %         %Get the dot product between the difference vector of each corner
        %         %of this voxel and the normal vector of the chosen nearest voxel
        %         %(the chosen voxel may not be the current voxel itself, as
        %         %neighbouring voxels are the same distance away)
        %         vox_dp = dot_products{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)));
        %         disp('Voxel corner dot products:');
        %         disp(num2str(vox_dp));
        %         disp(' ');
            end
        end %End i
        disp(['TOTAL number of octree cells with all control points having the same sign, at level ' num2str(lvl) ', before quantization: ' num2str(same_sign_cntr) '/' num2str(length(all_ctrlpts{lvl})/8) ' (' num2str((same_sign_cntr/(length(all_ctrlpts{lvl})/8))*100) '%)']);
        disp(['No. of cells with all 0 control points at level ' num2str(lvl) ': ' num2str(zero_cp_cntr)]);
        disp(' ');
        if (lvl == b + 1) && (same_sign_cntr > 0)
            %Plot voxels that have the same control point signs
            figure;
            %Original, input voxels
            scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, 'filled', 'MarkerFaceColor', 'b');
            hold on;
            %Voxels with same-sign control points
            scatter3(same_sign_voxels(:, 1), same_sign_voxels(:, 2), same_sign_voxels(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
            axis equal; axis off;
            title({'Voxels with Same-Sign Control Points at Encoder', 'Before Quantization'});
            legend('Original Voxels', 'Voxels with Same-Sign Control Points (BQ)', 'Location', 'best');
        end
    end %End lvl
end

% %------------- Visualization of Control Point Computation ----------------%
% 
% disp(' ');
% disp('Visualizing control point computation ...'); 
% disp('------------------------------------------------------------');
% 
% % Read in the input point cloud (assume PLY format)
% [~, ptcloud, ~] = plyRead(ptcloud_file);
% % 
% % %for lvl = 1:vis_levels_ctrlpts   
% for lvl = max_lvl   
%      %For each unique corner coordinate ...
%      %for c = 1:size(unique_coords{lvl}, 1)
%      for c = [36278 36285 36222 36186 36277 36284 36221 36182]
%         %Display the input point cloud
%         figure;
%         scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
%         axis equal; axis off;
%         hold on; 
%         %Display all the octree cells at this level
%         for j = 1:size(unique_coords{lvl}, 1)
%             %Get the list of vertex indices that share an edge with the 
%             %current vertex
%             edge_pts_indices = shared_edge_inds{lvl, j};
%             %Get the (x, y, z) coordinates corresponding to all the
%             %edge_pts_indices
%             edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%             %Connect the current vertex to all the edge_pts_coords with 
%             %straight lines
%             for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%                 plot3([unique_coords{lvl}(j, 1) edge_pts_coords(nv, 1)], [unique_coords{lvl}(j, 2) edge_pts_coords(nv, 2)], [unique_coords{lvl}(j, 3) edge_pts_coords(nv, 3)], 'Color', colours(lvl, :));
%                 hold on;
%             end
%         end
%         %Initialize a matrix to store all the corner coordinates of all the 
%         %octree cells at the current level, which share the current vertex 
%         shared_cell_coords = [];
%         %Get the indices of the shared octree cells (at the current level)
%         %for the current vertex
%         shared_cell_inds = shared_cells{lvl, c};  
%         %Get all the corner coordinates for all the shared cells
%         for sc = 1:length(shared_cell_inds)
%             if sc == 1
%                 shared_cell_coords = corner_coords{lvl}((shared_cell_inds(sc)*8 - 7):(shared_cell_inds(sc)*8), :);
%             else
%                 shared_cell_coords = cat(1, shared_cell_coords, corner_coords{lvl}((shared_cell_inds(sc)*8 - 7):(shared_cell_inds(sc)*8), :));
%             end
%         end
%         %Initialize a cell array to store the indices of the edge vertices
%         %for the shared cells of the current vertex
%         sc_edge_pts_indices = cell(size(shared_cell_coords, 1), 1);
%         %For each corner of the shared cells, figure out which other
%         %vertices in shared_cell_coords it is connected to
%         for scc = 1:size(shared_cell_coords, 1)
%             %Initialize a counter to keep track of how many shared corner
%             %coordinates have been found so far
%             sc_edge_pts_cntr = 1;
%             [x_same_rows, ~] = find(shared_cell_coords(scc, 1) == shared_cell_coords(:, 1));
%             [y_same_rows, ~] = find(shared_cell_coords(scc, 2) == shared_cell_coords(:, 2));
%             [z_same_rows, ~] = find(shared_cell_coords(scc, 3) == shared_cell_coords(:, 3));  
%             temp_cat = [x_same_rows; y_same_rows; z_same_rows];
%             %Only keep the vertex indices that appear twice in temp_cat
%             for v_ind = 1:length(temp_cat)
%                 if length(find(temp_cat == temp_cat(v_ind))) == 2
%                     sc_edge_pts_indices{scc}(sc_edge_pts_cntr) = temp_cat(v_ind);
%                     sc_edge_pts_cntr = sc_edge_pts_cntr + 1;
%                 end
%             end
%             %Keep only the unique indices in sc_edge_pts_indices{scc}
%             sc_edge_pts_indices{scc} = (unique(sc_edge_pts_indices{scc}))';
%         end
%         %Connect all the corner vertices of the shared octree cells, to
%         %other corner vertices on these cells, with which they share an 
%         %edge. The result will be a thick, black outline of the octree 
%         %cells at the current octree level, which share the current corner 
%         %vertex.
%         for c1 = 1:size(sc_edge_pts_indices, 1)
%             for c2 = 1:numel(sc_edge_pts_indices{c1})
%                 %Get the (x, y, z) coordinates of current two vertices that
%                 %will be connected
%                 xyz_c1 = shared_cell_coords(c1, :);
%                 xyz_c2 = shared_cell_coords(sc_edge_pts_indices{c1}(c2), :);
%                 h(2) = plot3([xyz_c1(1) xyz_c2(1)], [xyz_c1(2) xyz_c2(2)], [xyz_c1(3) xyz_c2(3)], 'k', 'LineWidth', 3);
%                 hold on;
%             end
%         end
%         title({['OCTREE LEVEL ' num2str(lvl) ':'], 'For Each Unique Cell Corner', 'Identifying Shared Octree Cells and', 'Computing the Signed Distance Value (Bezier Control Point)', 'from This Corner to the Nearest Occupied Voxel in the Shared Cells'}, 'Interpreter', 'none');
%         %Plot the current corner point on top of the input 3D point cloud 
%         %and the octree grid. If the control point associated with this
%         %corner is negative, plot the corner point in red with a black 
%         %outline; if the control point is positive, plot the control point
%         %in green with a black outline.
%         if (control_points{lvl}(c) < 0)
%             h(1) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
%             red_or_green = 'has -ive control point';
%         elseif (control_points{lvl}(c) > 0)
%             h(1) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
%             red_or_green = 'has +ive control point';
%         end
%         hold on;
%         %Colour in the nearest voxel (in blue, with a black outline) found
%         %for the current corner point. Note that the size of this voxel 
%         %will be exaggerated on this plot, in order to make it easier to 
%         %see.  
%         h(3) = scatter3(nearest_voxels{lvl}(c, 1), nearest_voxels{lvl}(c, 2), nearest_voxels{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
%         %text(nearest_voxels{lvl}(c, 1), nearest_voxels{lvl}(c, 2), nearest_voxels{lvl}(c, 3), ['(' num2str(nearest_voxels{lvl}(c, 1)) ', ' num2str(nearest_voxels{lvl}(c, 2)) ', ' num2str(nearest_voxels{lvl}(c, 3)) ')'], 'Color', 'b', 'FontSize', 18);
%         hold on;
%         %Also plot a vector to show the direction of the normal for the
%         %nearest voxel, and the difference vector between the nearest voxel
%         %and the corresponding corner point. Make the normal vector unit
%         %norm, then scale it by s, so that it is more easily visible. For
%         %the difference vector, just plot with its original size.
%         s = 80;
%         h(4) = quiver3(nearest_voxels{lvl}(c, 1), nearest_voxels{lvl}(c, 2), nearest_voxels{lvl}(c, 3), normal_nearest_vox{lvl}(c, 1)/norm(normal_nearest_vox{lvl}(c, :)), normal_nearest_vox{lvl}(c, 2)/norm(normal_nearest_vox{lvl}(c, :)), normal_nearest_vox{lvl}(c, 3)/norm(normal_nearest_vox{lvl}(c, :)), s, 'LineWidth', 2, 'Color', 'b', 'MaxHeadSize', 0.8);
%         hold on;
%         %h(5) = quiver3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), difference_vectors{lvl}(c, 1)/norm(difference_vectors{lvl}(c, :)), difference_vectors{lvl}(c, 2)/norm(difference_vectors{lvl}(c, :)), difference_vectors{lvl}(c, 3)/norm(difference_vectors{lvl}(c, :)), s, 'LineWidth', 2, 'Color', 'm', 'MaxHeadSize', 0.8);
%         h(5) = quiver3(nearest_voxels{lvl}(c, 1), nearest_voxels{lvl}(c, 2), nearest_voxels{lvl}(c, 3), difference_vectors{lvl}(c, 1), difference_vectors{lvl}(c, 2), difference_vectors{lvl}(c, 3), 'LineWidth', 2, 'Color', 'm', 'MaxHeadSize', 0.8);
%         legend(h, ['Current corner point (' red_or_green ')'], 'Shared octree cells (at current octree level)', 'Nearest voxel point', 'Normal direction for nearest voxel', 'Difference vector', 'Location', 'best');
%         hold off;
%      end %End current unique corner at level "lvl"
% end %End octree level "lvl"

% %---- Visualization of OT Cell Occupancy According to Control Points -----%
% 
% disp(' ');
% disp('Visualizing octree cell occupancy according to control points ...'); 
% disp('------------------------------------------------------------');
% 
% % Read in the input point cloud (assume PLY format)
% [~, ptcloud, ~] = plyRead(ptcloud_file);
% %For each octree level ...
% for lvl = [4]  
%     %Display the input point cloud
%     figure;
%     scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
%     axis equal; axis off;
%     hold on; 
%     %For each unique corner at this level ...
%     for c = 1:size(unique_coords{lvl}, 1)
%         %Display all the octree cells at this level ...
%         %Get the list of vertex indices that share an edge with the 
%         %current vertex
%         edge_pts_indices = shared_edge_inds{lvl, c};
%         %Get the (x, y, z) coordinates corresponding to all the
%         %edge_pts_indices
%         edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%         %Connect the current vertex to all the edge_pts_coords with 
%         %straight lines
%         for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%             plot3([unique_coords{lvl}(c, 1) edge_pts_coords(nv, 1)], [unique_coords{lvl}(c, 2) edge_pts_coords(nv, 2)], [unique_coords{lvl}(c, 3) edge_pts_coords(nv, 3)], 'Color', 'k');
%             hold on;
%         end
%         %Plot the current corner point on top of the input 3D point cloud 
%         %and the octree grid. If the control point associated with this
%         %corner is negative, plot the corner point in red with a black 
%         %outline; if the control point is positive, plot the control point
%         %in green with a black outline.
%         if (control_points{lvl}(c) < 0)
%             h_cp(1) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
%         elseif (control_points{lvl}(c) > 0)
%             h_cp(2) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
%         end
%         hold on;
%     end
%     hold off;
%     legend(h_cp, 'Has -ive control point', 'Has +ive control point', 'Location', 'best');
%     title(['Control Points at Octree Level ' num2str(lvl)]);
%     savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\ctrl_pts_visualization_lvl' num2str(lvl)]);
%     print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\ctrl_pts_visualization_lvl' num2str(lvl)], '-dpdf');
%     disp(['Level ' num2str(lvl) ' done']);
%     disp('------------------------------------------------------------');
% end

%------------------------- Wavelet Decomposition -------------------------%

disp(' ');
disp('----------------- Wavelet Decomposition --------------------');
disp(' ');

[wavelet_coeffs, reconstructed_control_points] = wavelet_analysis(debug_flag, myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, zero_threshold);
%[wavelet_coeffs, reconstructed_control_points] = wavelet_analysis_loop(myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, ptcloud_file);

if debug_flag == 1
    disp(' ');
    %For debugging purposes, check how many of the RECONSTRUCTED control 
    %points (out of the reconstructed control points for the UNIQUE 
    %occupied cell corners only) at each octree level, if any, are 0
    for lvl = start_lvl:1:max_lvl
        zero_ctrlpt_cnt = length(find(reconstructed_control_points{lvl} == 0));
        disp(['No. of 0 control points at level ' num2str(lvl) ', after quantization and reconstruction: ' num2str(zero_ctrlpt_cnt) '/' num2str(length(reconstructed_control_points{lvl}))]);
    end 
    disp(' ');

    %For debugging purposes, print out how many occupied octree cells at 
    %each octree level, if any, end up having all 8 of their RECONSTRUCTED 
    %control points (i.e., after quantization of control points at 
    %start_lvl and adding back quantized wavelet coefficients) with the 
    %same sign (+/-/0) ...

    %Accumulate the reconstructed control points for ALL the occupied 
    %octree cell corners at each level (not just the control points for the
    %unique corners) into one cell array
    all_ctrlpts = get_all_ctrlpts(reconstructed_control_points, ctrl_pts_pointers, start_lvl, max_lvl); 
    same_sign_voxels = [];
    disp(' ');
    for lvl = start_lvl:max_lvl
        same_sign_cntr = 0;
        cell_cntr = 0;
        zero_cp_cntr = 0;
        for i = 1:8:(length(all_ctrlpts{lvl}) - 7) 
            cell_cntr = cell_cntr + 1;
            %Get all 8 control points for the corners of the current cell
            current_ctrlpts = all_ctrlpts{lvl}(i:(i + 7));
            %Check if all control points of the current cell have the same 
            %sign, including the case where all the control points may be 0
            if (abs(sum(sign(current_ctrlpts))) == 8)||(~any(sign(current_ctrlpts)))
                same_sign_cntr = same_sign_cntr + 1;
                %disp(['Cell ' num2str(cell_cntr) ' has all control points with the same sign: ']);
                if ~any(sign(current_ctrlpts))
                    zero_cp_cntr = zero_cp_cntr + 1;
                end
                if lvl == b + 1
                    %Store this voxel, for debugging purposes
                    same_sign_voxels((size(same_sign_voxels, 1) + 1), 1:3) = ptcloud(cell_cntr, 1:3);
                end
                %Display the control points for all corners of this cell
                %disp(num2str(current_ctrlpts));
                %disp(' ');
                %Get the corner coordinates for each corner of this cell
                %cell_corners = corner_coords{lvl}((i:(i + lvl - 1)), 1:3);
                %disp('Cell corner coordinates:');
                %disp(num2str(cell_corners));
                %disp(' ');
            end
        end %End i
        disp(['TOTAL number of octree cells with all control points having the same sign, at level ' num2str(lvl) ', after quantization and reconstruction: ' num2str(same_sign_cntr) '/' num2str(length(all_ctrlpts{lvl})/8) ' (' num2str((same_sign_cntr/(length(all_ctrlpts{lvl})/8))*100) '%)']);
        disp(['No. of cells with all 0 control points at level ' num2str(lvl) ': ' num2str(zero_cp_cntr)]);
        disp(' ');
        if (lvl == b + 1) && (same_sign_cntr > 0)
            %Plot voxels that have the same control point signs
            figure;
            %Original, input voxels
            scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, 'filled', 'MarkerFaceColor', 'b');
            hold on;
            %Voxels with same-sign control points
            scatter3(same_sign_voxels(:, 1), same_sign_voxels(:, 2), same_sign_voxels(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
            axis equal; axis off;
            title({'Voxels with Same-Sign Control Points at Encoder', 'After Quantization and Reconstruction'});
            legend('Original Voxels', 'Voxels with Same-Sign Control Points (AQ)', 'Location', 'best');
        end
    end %End lvl
end

%----------------------------- Distributions -----------------------------%

% if debug_flag == 1
%     %Plot the distribution of reconstructed control points (low-pass
%     %coefficients) at different octree levels
%     for lvl = start_lvl:max_lvl
%         figure;
%         histogram(reconstructed_control_points{lvl});
%         title(['Histogram of (Encoder-Reconstructed) Control Points at Octree Level ' num2str(lvl)]);
%         %Save the above histogram as a MATLAB figure and as a PDF image in 
%         %our network directory (NB: The '-bestfit' option maximizes the 
%         %size of the figure to fill the page, but preserves the aspect 
%         %ratio of the figure. The figure might not fill the entire page. 
%         %This option leaves a minimum page margin of .25 inches).
%         savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_ctrlpts_lvl' num2str(lvl)]);
%         print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_ctrlpts_lvl' num2str(lvl)], '-dpdf');
%     end
%     %Plot the distribution of quantized wavelet coefficients at different
%     %octree levels
%     for lvl = start_lvl:(max_lvl - 1)
%         figure;
%         histogram(wavelet_coeffs{lvl + 1});
%         title(['Histogram of Quantized Wavelet Coefficients between Octree Levels ' num2str(lvl) ' and ' num2str(lvl + 1)]);
%         %Save the above histogram as a MATLAB figure and as a PDF image in 
%         %our network directory (NB: The '-bestfit' option maximizes the 
%         %size of the figure to fill the page, but preserves the aspect 
%         %ratio of the figure. The figure might not fill the entire page. 
%         %This option leaves a minimum page margin of .25 inches).
%         savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_wavcfs_lvls' num2str(lvl) '-' num2str(lvl + 1)]);
%         print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_wavcfs_lvls' num2str(lvl) '-' num2str(lvl + 1)], '-dpdf');
%     end
% end

%---------------- Checking for Zero Wavelet Coefficients -----------------%

disp(' ');
disp('------- Checking for All Zero Wavelet Coefficients  --------');
disp(' ');

%Cell array to store the indices of the occupied octree cells at each level
%that have zero wavelet coefficients (quantized symbols) on all of their 
%corners
all_zero_wav_cfs = cell(size(wavelet_coeffs));
%Cell array to store the midpoints of the occupied octree cells at each
%level that have zero wavelet coefficients on all of their corners (only
%needed for display purposes, later)
%all_zero_cell_midpoints = cell(size(wavelet_coeffs));
%Cell array to store the indices of ctrl_pts_pointers for the occupied
%octree cells at each level that have zero wavelet coefficients on all of
%their corners (only needed for display purposes, later)
%all_zero_cell_ctrlpts_ptrs = cell(size(wavelet_coeffs));

%At each octree level, check which occupied octree cells (if any) have zero
%wavelet coefficients on all of their corners
for lvl = (start_lvl + 1):max_lvl
    if debug_flag == 1
        disp(['Processing octree level ' num2str(lvl) ' ...']);
    end
    %Counter for number of occupied cells at this level, which contain all
    %zero wavelet coefficients
    zw_cntr = 1;
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Get the wavelet coefficient for each of the 8 corners of this cell
        current_wavelet_coeffs = wavelet_coeffs{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %If the quantized wavelet coefficients at all the corners of this 
        %cell are 0 ...
        if isempty(find((current_wavelet_coeffs ~= 0), 1))
        %If the wavelet coefficients at all the corners of this cell are
        %near 0 (i.e., fit within +/- zero_threshold of 0) ...
%         if isempty(find((abs(current_wavelet_coeffs) > zero_threshold), 1))
            %disp(['All zero wavelet coefficients for occupied cell ' num2str(occ_cell)]);
            all_zero_wav_cfs{lvl}(zw_cntr) = occ_cell; 
            %Get the midpoints of the current cell (only needed for display
            %purposes, later)
            %all_zero_cell_midpoints{lvl}(zw_cntr, 1:3) = mean(corner_coords{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :), 1);
            %Store the control points pointers for each corner of this cell
            %(only needed for display purposes, later)
            %all_zero_cell_ctrlpts_ptrs{lvl}(zw_cntr, 1:8) = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8));
            zw_cntr = zw_cntr + 1;
        end
    end
    if debug_flag == 1
        disp(['TOTAL no. of occupied octree cells with all zero wavelet coefficients at this level: ' num2str(length(all_zero_wav_cfs{lvl}))]);
        disp('------------------------------------------------------------');
    end
end

% %---------- Visualization of Zero Wavelet Coefficient Locations ----------%
% 
% disp(' ');
% disp('Visualizing locations of zero wavelet coefficients on octree cell corners ...'); 
% disp('------------------------------------------------------------');
% 
% % Read in the input point cloud (assume PLY format)
% [~, ptcloud, ~] = plyRead(ptcloud_file);
% %For each octree level ...
% for lvl = [4, 5]  
%     %Display the input point cloud
%     figure;
%     scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
%     axis equal; axis off;
%     hold on; 
%     %For each unique corner at this level ...
%     for c = 1:size(unique_coords{lvl}, 1)
%         %Display all the octree cells at this level ...
%         %Get the list of vertex indices that share an edge with the 
%         %current vertex
%         edge_pts_indices = shared_edge_inds{lvl, c};
%         %Get the (x, y, z) coordinates corresponding to all the
%         %edge_pts_indices
%         edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%         %Connect the current vertex to all the edge_pts_coords with 
%         %straight lines
%         for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%             plot3([unique_coords{lvl}(c, 1) edge_pts_coords(nv, 1)], [unique_coords{lvl}(c, 2) edge_pts_coords(nv, 2)], [unique_coords{lvl}(c, 3) edge_pts_coords(nv, 3)], 'Color', 'k');
%             hold on;
%         end
%         %If this corner has a zero wavelet coefficient, plot the corner 
%         %point in red with a black outline; otherwise, don't plot the
%         %corner point
%         if (wavelet_coeffs{lvl}(c) == 0)
%             h_zw(1) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 50, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
%         end
%         hold on;
%     end
%     if ~isempty(all_zero_wav_cfs{lvl})
%         %For each octree cell at this level, which has all zero wavelet
%         %coefficients
%         for ocz = 1:length(all_zero_wav_cfs{lvl})
%             %Get the coordinates of all 8 corners of this cell
%             current_cell_coords = corner_coords{lvl}(((all_zero_wav_cfs{lvl}(ocz)*8 - 7):(all_zero_wav_cfs{lvl}(ocz)*8)), :);
%             %Connect the corners that share an edge, with thick red lines
%             %(see p86 in Logbook1 for the corner ordering used below)
%             %Connect corner 1 to corners 2, 4, and 5
%             for cnr2 = [2, 4, 5]
%                 h_zw(2) = plot3([current_cell_coords(1, 1) current_cell_coords(cnr2, 1)], [current_cell_coords(1, 2) current_cell_coords(cnr2, 2)], [current_cell_coords(1, 3) current_cell_coords(cnr2, 3)], 'Color', 'r', 'LineWidth', 3);
%                 hold on;
%             end
%             %Connect corner 2 to corners 3 and 6
%             for cnr2 = [3, 6]
%                 h_zw(2) = plot3([current_cell_coords(2, 1) current_cell_coords(cnr2, 1)], [current_cell_coords(2, 2) current_cell_coords(cnr2, 2)], [current_cell_coords(2, 3) current_cell_coords(cnr2, 3)], 'Color', 'r', 'LineWidth', 3);
%                 hold on;
%             end
%             %Connect corner 3 to corners 4 and 7
%             for cnr2 = [4, 7]
%                 h_zw(2) = plot3([current_cell_coords(3, 1) current_cell_coords(cnr2, 1)], [current_cell_coords(3, 2) current_cell_coords(cnr2, 2)], [current_cell_coords(3, 3) current_cell_coords(cnr2, 3)], 'Color', 'r', 'LineWidth', 3);
%                 hold on;
%             end
%             %Connect corner 4 to corner 8
%             h_zw(2) = plot3([current_cell_coords(4, 1) current_cell_coords(8, 1)], [current_cell_coords(4, 2) current_cell_coords(8, 2)], [current_cell_coords(4, 3) current_cell_coords(8, 3)], 'Color', 'r', 'LineWidth', 3);
%             hold on;
%             %Connect corner 5 to corners 6 and 8
%             for cnr2 = [6, 8]
%                 h_zw(2) = plot3([current_cell_coords(5, 1) current_cell_coords(cnr2, 1)], [current_cell_coords(5, 2) current_cell_coords(cnr2, 2)], [current_cell_coords(5, 3) current_cell_coords(cnr2, 3)], 'Color', 'r', 'LineWidth', 3);
%                 hold on;
%             end
%             %Connect corner 6 to corner 7
%             h_zw(2) = plot3([current_cell_coords(6, 1) current_cell_coords(7, 1)], [current_cell_coords(6, 2) current_cell_coords(7, 2)], [current_cell_coords(6, 3) current_cell_coords(7, 3)], 'Color', 'r', 'LineWidth', 3);
%             hold on;
%             %Connect corner 7 to corner 8
%             h_zw(2) = plot3([current_cell_coords(7, 1) current_cell_coords(8, 1)], [current_cell_coords(7, 2) current_cell_coords(8, 2)], [current_cell_coords(7, 3) current_cell_coords(8, 3)], 'Color', 'r', 'LineWidth', 3);
%             hold on;
%             %Display the occupied cell index in the centre of each occupied 
%             %cell at this level that contains all zero wavelet coefficients
%             text(all_zero_cell_midpoints{lvl}(ocz, 1), all_zero_cell_midpoints{lvl}(ocz, 2), all_zero_cell_midpoints{lvl}(ocz, 3), num2str(all_zero_wav_cfs{lvl}(ocz)), 'Color', 'b', 'FontSize', 18);
%             hold on;
%         end
%     end   
%     hold off;
%     if ~isempty(all_zero_wav_cfs{lvl})
%         legend(h_zw, 'Has zero wavelet coefficient', 'Cell with all zero wavelet coefficients', 'Location', 'best');
%     else
%         legend(h_zw, 'Has zero wavelet coefficient', 'Location', 'best');
%     end
%     title(['Locations of Zero Wavelet Coefficients at Octree Level ' num2str(lvl)]);
%     savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\zero_wavelet_coeffs_lvl' num2str(lvl)]);
%     print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\zero_wavelet_coeffs_lvl' num2str(lvl)], '-dpdf');
%     disp(['Level ' num2str(lvl) ' done']);
%     disp('------------------------------------------------------------');
% end

%------------------------------- Pruning  --------------------------------%

%Prune octree and wavelet coefficient tree, if pruning is required
if prune_flag == 1
    disp(' ');
    disp('-------------------- Pruning Octree ------------------------');
    disp(' ');

    [pruned_occupancy_codes, post_pruning_array, toprune, toprune2] = prune_octree(debug_flag, myOT, all_zero_wav_cfs, start_lvl, max_lvl, b, reconstructed_control_points, ctrl_pts_pointers);

    disp(' ');
    disp('------------ Pruning Wavelet Coefficient Tree --------------');
    disp(' ');

    pruned_wavelet_coeffs = prune_wavelet_coeff_tree(debug_flag, wavelet_coeffs, toprune, toprune2, ctrl_pts_pointers, myOT, start_lvl, max_lvl, b);
end

% %The below is for debugging only: prune the reconstructed control points to
% %match the control points that will be reconstructed at the decoder, to
% %make sure that they are exactly the same
% pruned_reconstructed_control_points = cell(size(reconstructed_control_points));
% for i = start_lvl:(max_lvl - 1)
%     if isempty(toprune2{i})
%         %The below will represent only the control points for the UNIQUE
%         %corners
%         pruned_reconstructed_control_points{i} = reconstructed_control_points{i};
%     else
%         first_inds = toprune2{i}.*8 - 7;
%         last_inds = toprune2{i}.*8;
%         all_inds = [];
%         for j = 1:length(first_inds)
%             all_inds((end + 1):(end + 8), 1) = first_inds(j):last_inds(j);
%         end
%         pruned_reconstructed_control_points{i} = reconstructed_control_points{i}(ctrl_pts_pointers{i});
%         %The below will represent ALL the control points at each level
%         %after pruning (except the voxel level), not just the unique corner
%         %control points. NOTE that this differs from the control points
%         %stored in pruned_reconstructed_control_points at levels where
%         %toprune2 is empty, above, which are stored only for the unique
%         %corners. To compare the below with control points reconstructed
%         %at the decoder, expand the decoder-reconstructed control points by
%         %using ctrl_pts_pointers, to get the control points for ALL corners
%         %at the corresponding octree level.
%         pruned_reconstructed_control_points{i}(all_inds) = [];
%     end
% end
% save('pruned_reconstructed_control_points', 'pruned_reconstructed_control_points');
 
%------------------------------- Encoding --------------------------------%

disp(' ');
disp('--------------- Encoding for Transmission ------------------');
disp(' ');

start_compute_bitrates_time = tic;

%---- Occupancy Codes ----

%Get the octree occupancy codes that will be transmitted to the decoder
occupancy_codes_forDec = cell(b, 1);
for i = 1:(max_lvl - 1)
    if prune_flag == 1
        occupancy_codes_forDec{i} = pruned_occupancy_codes{i};
    else
        occupancy_codes_forDec{i} = myOT.OccupancyCode{i};
    end
end
%Concatenate all of the occupancy codes (decimal values) at all octree 
%levels from the root to (max_lvl - 1), into one long array
oc_cntr = 1;
occ_codes_array = [];
for l = 1:(max_lvl - 1)
    occ_codes_array(oc_cntr:(oc_cntr + numel(occupancy_codes_forDec{l}) - 1)) = occupancy_codes_forDec{l};
    oc_cntr = oc_cntr + numel(occupancy_codes_forDec{l});
end
if debug_flag == 1
    %Plot a histogram of the occupancy codes inside occ_codes_array
    figure;
    histogram(occ_codes_array);
    if prune_flag == 1
        title(['Histogram of Pruned Octree Occupancy Codes from Level 1-' num2str(max_lvl - 1)]);
    elseif prune_flag == 0
        title(['Histogram of Octree Occupancy Codes from Level 1-' num2str(max_lvl - 1)]);
    end
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_occ_codes_lvl1-' num2str(max_lvl - 1)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_occ_codes_lvl1-' num2str(max_lvl - 1)], '-dpdf');
end
%Compute the number of bits required for these occupancy codes
bits_occ_codes_persymbol = entropy_calc(occ_codes_array);    %Avg. minimum no. of bits per symbol
bits_occ_codes = bits_occ_codes_persymbol*length(occ_codes_array);    %Total no. of bits for all symbols
disp(['Total entropy bits for occupancy codes: ' num2str(bits_occ_codes) ' (' num2str(bits_occ_codes_persymbol) ' bits per symbol)']);
disp(['Occupancy codes bpv (bits per input voxel): ' num2str(bits_occ_codes/size(ptcloud, 1))]);
disp(' ');

if debug_flag == 1
    %Do the below for comparison only, of the number of bits required for
    %unpruned vs pruned occupancy codes
    if prune_flag == 1
        %Get the unpruned occupancy codes that would be transmitted to the 
        %decoder
        occupancy_codes_wout_pruning = cell(b, 1);
        for i = 1:(max_lvl - 1)
            occupancy_codes_wout_pruning{i} = myOT.OccupancyCode{i};
        end
        %Concatenate all of the occupancy codes (decimal values) at all
        %octree levels from the root to (max_lvl - 1), into one long array
        oc_cntr_wout_pruning = 1;
        occ_codes_wout_pruning = [];
        for l = 1:(max_lvl - 1)
            occ_codes_wout_pruning(oc_cntr_wout_pruning:(oc_cntr_wout_pruning + numel(occupancy_codes_wout_pruning{l}) - 1)) = occupancy_codes_wout_pruning{l};
            oc_cntr_wout_pruning = oc_cntr_wout_pruning + numel(occupancy_codes_wout_pruning{l});
        end
        %Compute the number of bits that would be required for the unpruned
        %occupancy codes
        disp('******** FOR COMPARISON ONLY ********');
        disp(' ');
        bits_unpruned_occ_codes_persymbol = entropy_calc(occ_codes_wout_pruning);    %Avg. minimum no. of bits per symbol
        bits_unpruned_occ_codes = bits_unpruned_occ_codes_persymbol*length(occ_codes_wout_pruning);    %Total no. of bits for all symbols
        disp(['Total entropy bits for UNpruned occupancy codes: ' num2str(bits_unpruned_occ_codes) ' (' num2str(bits_unpruned_occ_codes_persymbol) ' bits per symbol)']);
        disp(['UNpruned occupancy codes bpv (bits per input voxel): ' num2str(bits_unpruned_occ_codes/size(ptcloud, 1))]);
        disp(' ');
        disp('*************************************');
        disp(' ');
    end
end

%---- Post-Pruning Array (only if pruning has been done) ----

if prune_flag == 1    
    %Get the (parts of the) post_pruning_array that will be transmitted to 
    %the decoder
    post_pruning_array_forDec = cell(b, 1);
    if max_lvl == b + 1
        for i = 1:(max_lvl - 1)
            post_pruning_array_forDec{i} = post_pruning_array{i};
        end
    elseif max_lvl < b + 1
        for i = 1:max_lvl
            post_pruning_array_forDec{i} = post_pruning_array{i};
        end
    end
    %Concatenate all of the bits from post_pruning_array_forDec, at all 
    %octree levels from the root to (max_lvl - 1), into one long array
    pp_cntr = 1;
    pp_array = [];
    if max_lvl == b + 1
        for l = 1:(max_lvl - 1)
            pp_array(pp_cntr:(pp_cntr + numel(post_pruning_array_forDec{l}) - 1)) = post_pruning_array_forDec{l};
            pp_cntr = pp_cntr + numel(post_pruning_array_forDec{l});
        end
    elseif max_lvl < b + 1
        for l = 1:max_lvl
            pp_array(pp_cntr:(pp_cntr + numel(post_pruning_array_forDec{l}) - 1)) = post_pruning_array_forDec{l};
            pp_cntr = pp_cntr + numel(post_pruning_array_forDec{l});
        end
    end
    if debug_flag == 1
        %Plot a histogram of the bits inside pp_array
        figure;
        histogram(pp_array);
        title({'Histogram of Post-Pruning Array Bits', 'Indicating Leaf/Non-Leaf Octree Cells', ['from Level 1-' num2str(max_lvl - 1)]});
        %Save the above histogram as a MATLAB figure and as a PDF image in 
        %our network directory (NB: The '-bestfit' option maximizes the 
        %size of the figure to fill the page, but preserves the aspect 
        %ratio of the figure. The figure might not fill the entire page. 
        %This option leaves a minimum page margin of .25 inches).
        savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_pp_bits_lvl1-' num2str(max_lvl - 1)]);
        print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_pp_bits_lvl1-' num2str(max_lvl - 1)], '-dpdf');
    end
    %Compute the number of bits required for pp_array
    bits_pp_array_persymbol = entropy_calc(pp_array);    %Avg. minimum no. of bits per symbol
    bits_pp_array = bits_pp_array_persymbol*length(pp_array);    %Total no. of bits for all symbols
    disp(['Total entropy bits for post-pruning array: ' num2str(bits_pp_array) ' (' num2str(bits_pp_array_persymbol) ' bits per symbol)']);
    disp(['pp_array bpv (bits per input voxel): ' num2str(bits_pp_array/size(ptcloud, 1))]);
    disp(' ');
    
    %Set the output variable
    varargout{1} = post_pruning_array_forDec;    
else
    varargout = {};
end %End check if prune_flag == 1 

%---- Control Points at start_lvl ----

%Get the quantized control points at start_lvl, which will be transmitted 
%to the decoder
%rec_ctrlpts_forDec = reconstructed_control_points{start_lvl};
rec_ctrlpts_forDec = quantize_uniform_scalar(control_points{start_lvl}, q_stepsize);
if debug_flag == 1
    %Plot a histogram of the quantized control points inside 
    %rec_ctrlpts_forDec
    figure;
    histogram(rec_ctrlpts_forDec);
    title(['Histogram of Quantized Control Points at Start Level (Level ' num2str(start_lvl) ')']);
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_ctrlpts_lvl' num2str(start_lvl)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_ctrlpts_lvl' num2str(start_lvl)], '-dpdf');
end
%Compute the number of bits required for these control points
bits_ctrlpts_persymbol = entropy_calc(rec_ctrlpts_forDec);  %Avg. minimum no. of bits per symbol
bits_ctrlpts = bits_ctrlpts_persymbol*length(rec_ctrlpts_forDec);    %Total no. of bits for all symbols
disp(['Total entropy bits for control points at start_lvl: ' num2str(bits_ctrlpts) ' (' num2str(bits_ctrlpts_persymbol) ' bits per symbol)']);
disp(['start_lvl control points bpv (bits per input voxel): ' num2str(bits_ctrlpts/size(ptcloud, 1))]);
disp(' ');

%---- Wavelet Coefficients ----

%Get the quantized wavelet coefficients that will be transmitted to the
%decoder
wavelet_coeffs_forDec = cell(b, 1);
for j = (start_lvl + 1):max_lvl
    if prune_flag == 1
        wavelet_coeffs_forDec{j} = pruned_wavelet_coeffs{j}; 
    else
        wavelet_coeffs_forDec{j} = wavelet_coeffs{j};
    end
end
%Concatenate all of the wavelet coefficients at all octree levels from 
%(start_lvl + 1) to max_lvl, into one long array
wcf_cntr = 1;
wavelet_cfs_array = [];
for l = (start_lvl + 1):max_lvl
    wavelet_cfs_array(wcf_cntr:(wcf_cntr + numel(wavelet_coeffs_forDec{l}) - 1)) = wavelet_coeffs_forDec{l};
    wcf_cntr = wcf_cntr + numel(wavelet_coeffs_forDec{l});
end
if debug_flag == 1
    %Plot a histogram of the quantized wavelet coefficients inside 
    %wavelet_cfs_array
    figure;
    histogram(wavelet_cfs_array);
    if prune_flag == 1
        title({'Histogram of Quantized Pruned Wavelet Coefficients', ['from Level ' num2str((start_lvl + 1)) '-' num2str(max_lvl)]});
    elseif prune_flag == 0
        title(['Histogram of Quantized Wavelet Coefficients from Level ' num2str((start_lvl + 1)) '-' num2str(max_lvl)]);
    end
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_wavcfs_lvl' num2str(start_lvl + 1) '-' num2str(max_lvl)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_wavcfs_lvl' num2str(start_lvl + 1) '-' num2str(max_lvl)], '-dpdf');
end
%Compute the number of bits required for these wavelet coefficients
bits_wavelet_cfs_persymbol = entropy_calc(wavelet_cfs_array);  %Avg. minimum no. of bits per symbol
bits_wavelet_cfs = bits_wavelet_cfs_persymbol*length(wavelet_cfs_array);    %Total no. of bits for all symbols
disp(['Total entropy bits for wavelet coefficients: ' num2str(bits_wavelet_cfs) ' (' num2str(bits_wavelet_cfs_persymbol) ' bits per symbol)']);
disp(['Wavelet coefficients bpv (bits per input voxel): ' num2str(bits_wavelet_cfs/size(ptcloud, 1))]);
disp(' ');

if debug_flag == 1
    %Do the below for comparison only, of the number of bits required for
    %unpruned vs pruned wavelet coefficients
    if prune_flag == 1
        %Get the unpruned quantized wavelet coefficients that would be 
        %transmitted to the decoder
        wavelet_coeffs_wout_pruning = cell(b, 1);
        for j = (start_lvl + 1):max_lvl
            wavelet_coeffs_wout_pruning{j} = wavelet_coeffs{j};
        end
        %Concatenate all of the wavelet coefficients at all octree levels
        %from (start_lvl + 1) to max_lvl, into one long array
        wcf_cntr_wout_pruning = 1;
        wcf_array_wout_pruning = [];
        for l = (start_lvl + 1):max_lvl
            wcf_array_wout_pruning(wcf_cntr_wout_pruning:(wcf_cntr_wout_pruning + numel(wavelet_coeffs_wout_pruning{l}) - 1)) = wavelet_coeffs_wout_pruning{l};
            wcf_cntr_wout_pruning = wcf_cntr_wout_pruning + numel(wavelet_coeffs_wout_pruning{l});
        end
        %For comparison only, compute the number of bits that would be 
        %required for the unpruned wavelet coefficients
        disp('******** FOR COMPARISON ONLY ********');
        disp(' ');
        bits_unpruned_wcf_persymbol = entropy_calc(wcf_array_wout_pruning);    %Avg. minimum no. of bits per symbol
        bits_unpruned_wcf = bits_unpruned_wcf_persymbol*length(wcf_array_wout_pruning);    %Total no. of bits for all symbols
        disp(['Total entropy bits for UNpruned wavelet coefficients: ' num2str(bits_unpruned_wcf) ' (' num2str(bits_unpruned_wcf_persymbol) ' bits per symbol)']);
        disp(['UNpruned wavelet coefficients bpv (bits per input voxel): ' num2str(bits_unpruned_wcf/size(ptcloud, 1))]);
        disp(' ');
        disp('*************************************');
        disp(' ');
    end
end

%---- TOTALS (Geometry) ----

%Compute total geometry bits for transmission
if prune_flag == 1
    total_geom_bits = bits_occ_codes + bits_pp_array + bits_ctrlpts + bits_wavelet_cfs;
else
    total_geom_bits = bits_occ_codes + bits_ctrlpts + bits_wavelet_cfs;
end
total_geom_bpv = total_geom_bits/size(ptcloud, 1);
disp(['TOTAL bits: ' num2str(total_geom_bits)]);
disp(['TOTAL bpv: ' num2str(total_geom_bpv)]);

compute_bitrates_time = toc(start_compute_bitrates_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all bitrates for transmission: ' num2str(compute_bitrates_time) ' seconds']);
disp('************************************************************');

total_encoder_time = toc(start_enc_time);

disp(' ');
disp('------------------- ENCODER FINISHED -----------------------');
disp(' ');

disp(['TOTAL ENCODER TIME: ' num2str(total_encoder_time) ' seconds']);
disp(' ');

% %For debugging only: check which, if any, of the voxels that were not
% %reconstructed at the decoder, have control points with all the same signs
% load('test_vox_diffs.mat');
% same_sign_cntr = 0;
% for i = 1:size(test_vox_diffs, 1)
%     diff_temp = test_vox_diffs(i, :) - ptcloud(:, 1:3);
%     vox_ind = find(sum(abs(diff_temp), 2) == 0);
%     curr_ctrlpts = all_ctrlpts{8}((8*vox_ind - 7):(8*vox_ind)); %All control points, before pruning
%     if (sum(sign(curr_ctrlpts)) == length(curr_ctrlpts))||(sum(sign(curr_ctrlpts)) == -length(curr_ctrlpts))
%         same_sign_cntr = same_sign_cntr + 1;
%         disp(['Vox. ' num2str(i) ' (' num2str(test_vox_diffs(i, :)) ') has all control points with the same sign: ']);
%         disp(num2str(curr_ctrlpts));
%     end
% end
% disp(['Number of not-reconstructed voxels with all control points having the same sign: ' num2str(same_sign_cntr) '/' num2str(size(test_vox_diffs, 1))]);
% disp(' ');










