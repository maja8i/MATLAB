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

function [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points, varargout] = Bezier_volumes_encoder(ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag)

disp(' ');
disp('============================================================');
disp('                   ENCODER RUNNING ...');
disp('============================================================');
disp(' ');

%Initialize a cell array to store the corner coordinates of each occupied
%cell at each level of the octree
corner_coords = cell((b + 1), 1);
%Initialize a cell array to store only the unique corner coordinates at
%each octree level
unique_coords = cell((b + 1), 1);
%Initialize a cell array to store the Bezier control points (signed
%distances) associated with the unique corner coordinates at each octree 
%level
control_points = cell((b + 1), 1);
%Initialize a cell array of "pointers": there will be 8 pointers per
%occupied octree cell (1 per corner) at each octree level, which will point
%to the Bezier control point in control_points, associated with that
%corner. So, for each octree level "lvl", there will be 
%myOT.NodeCount(lvl) pointer arrays containining 8 elements each.  
ctrl_pts_pointers = cell((b + 1), 1);

%-------------------------- Octree Construction --------------------------%

disp('-------------------- Octree Construction -------------------');
disp(' ');

%Construct an octree of depth (b + 1), for the input point cloud
%[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted, original_points_per_voxel, original_normals_per_voxel] = construct_octree(ptcloud_file, b, ptcloud_name);
[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted] = construct_octree(ptcloud_file, b);

disp(' ');
disp('-------------- Extracting Voxel Coordinates ----------------');
disp(' ');

%Extract the set of occupied voxel coordinates (x, y, z) at all levels of
%the octree myOT 
%[~, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_centroids] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted, original_points_per_voxel, original_normals_per_voxel);
[~, occupied_voxel_coords, occupied_voxel_normals] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, normals_sorted);

disp(' ');
disp('-------------- Computing Corner Coordinates ----------------');
disp(' ');

%For each octree level ...
tic;
%for lvl = 1:(b + 1) 
for lvl = 1:max_lvl 

    disp(['Processing octree level ' num2str(lvl) ':']); 
    
    disp('Computing corner coordinates for each occupied cell ...');
    disp('------------------------------------------------------------');
    
    %Find the (x, y, z) coordinates of the origin of each occupied octree
    %cell at the current level (origin is at the bottom left-hand corner 
    %farthest from the viewer)
    corners1 = double(myOT.SpatialIndex{lvl})*(2^(b + 1 - lvl)) - [0.5 0.5 0.5];
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
cornercoords_time = toc;
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
%Initialize a cell array to store, for each unique corner vertex at each
%octree level, a list of indices of the corner vertices that share an edge
%with that vertex. 
shared_edge_inds = cell(b, 1);
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
for lvl = [4, 5]
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
    for v = 1:size(unique_coords{lvl}, 1)
        %Initialize an array to store the indices of the shared vertices,
        %for the current vertex
        edge_pts_indices = [];
        %Initialize a counter to keep track of how many vertex indices
        %corresponding to shared edges have been found
        edge_pts_indices_cntr = 1;
        %Get the coordinates of the current vertex
        current_vtx = unique_coords{lvl}(v, :);
        %Find all the vertices (amongst the list of unique corner vertices
        %at the current octree level) that share an edge with current_vtx
        [x_same_rows, ~] = find(current_vtx(1) == unique_coords{lvl}(:, 1));
        x_same_rows(find(x_same_rows == v)) = [];   %Remove the current vertex index
        [y_same_rows, ~] = find(current_vtx(2) == unique_coords{lvl}(:, 2));
        y_same_rows(find(y_same_rows == v)) = [];   %Remove the current vertex index
        [z_same_rows, ~] = find(current_vtx(3) == unique_coords{lvl}(:, 3));
        z_same_rows(find(z_same_rows == v)) = [];   %Remove the current vertex index
        temp_cat = [x_same_rows; y_same_rows; z_same_rows];
        for v_ind = 1:length(temp_cat)
            %Only keep the vertex indices that appear twice in temp_cat
            if length(find(temp_cat == temp_cat(v_ind))) == 2
                edge_pts_indices(edge_pts_indices_cntr) = temp_cat(v_ind);
                edge_pts_indices_cntr = edge_pts_indices_cntr + 1;
            end
        end
        %Keep only the unique indices in edge_pts_indices
        edge_pts_indices = (unique(edge_pts_indices))';
        %Get the (x, y, z) coordinates corresponding to all the
        %edge_pts_indices
        edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%         %Connect the current_vtx to all the edge_pts_coords with straight
%         %lines
%         for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%             h_ot(lvl) = plot3([current_vtx(1) edge_pts_coords(nv, 1)], [current_vtx(2) edge_pts_coords(nv, 2)], [current_vtx(3) edge_pts_coords(nv, 3)], 'Color', colours(lvl, :));
%             hold on;
%         end
        %Store the current list of edge indices inside shared_edge_inds,
        %for future reference
        shared_edge_inds{lvl, v} = edge_pts_indices;
    end
end
% hold off;
% legend(h_ot, otplot_legend, 'Location', 'best');
% title({'Octree Subdivision and Computing Corner Coordinates of Octree Cells', '(ENCODER)'}, 'Interpreter', 'none');

%-------------------- Computing Bezier Control Points --------------------%

disp(' ');
disp('------------ Computing Bezier Control Points ---------------');
disp(' ');

%Choose the highest octree level at which you wish the control points to be
%computed
%max_lvl = 8;

%Initialize cell array to store the indices of octree cells that share each
%unique corner, at each octree level
shared_cells = cell((b + 1), 1);
%Initialize a cell array to store the (x, y, z) coordinates of the nearest
%voxel found, for each unique corner coordinate, at each octree level
nearest_voxels = cell((b + 1), 1);
%Initialize a cell array to store the (x, y, z) coordinates of the normal
%vector of the nearest voxel for each unique corner coordinate, at each
%octree level
normal_nearest_vox = cell((b + 1), 1);
%Initialize a cell array to store the difference vector between the normal
%vector of the nearest voxel for each unique corner coordinate, and that 
%corner coordinate itself, at each octree level
difference_vectors = cell((b + 1), 1);
%Initialize a cell array to store the dot products between the difference
%vectors and normal vectors at each octree level
dot_products = cell((b + 1), 1);
%Initialize a cell array to store the absolute values of the minimum
%Euclidean distances (unsigned control points) for each unique corner at
%each octree level
min_euclid_dist = cell((b + 1), 1);

tic;
%for lvl = 1:(b + 1)
for lvl = 1:max_lvl
    disp(['Computing Bezier control points for octree level ' num2str(lvl) ':']);
    %Initialize a cell array to store the indices of any equidistant 
    %voxels found for a given corner at the voxel level
    if lvl == b + 1
        equidistant_voxel_inds = cell(size(unique_coords{b + 1}, 1), 1); 
    end
    %For each unique corner coordinate, find the nearest occupied voxel in
    %any of the occupied octree cells at the current level that share this 
    %corner. Do this by measuring the Euclidean distance from the corner to 
    %the (centres of) voxels.
    disp('Computing signed distance to nearest occupied voxel, for each unique corner coordinate ...');  
    for c = 1:size(unique_coords{lvl}, 1)
        %Find out which occupied octree cells share this corner
        shared_cells{lvl, c} = ceil(find((c - ctrl_pts_pointers{lvl}) == 0)./8);    %Faster than using "ismember"                
        %Get the (x, y, z) coordinates of all the occupied voxels in the 
        %cells found above, and the normals for each of these voxels
        current_occupied_voxel_coords = occupied_voxel_coords{lvl, shared_cells{lvl, c}(1)};
        current_occupied_voxel_normals = occupied_voxel_normals{lvl, shared_cells{lvl, c}(1)};
        if length(shared_cells{lvl, c}) > 1
            %Concatenate the position and normal coordinates from the cells 
            %that share the current corner (concatenate along rows, so that
            %we maintain the 3-column structure)
            current_occupied_voxel_coords = cat(1, current_occupied_voxel_coords, occupied_voxel_coords{lvl, shared_cells{lvl, c}(2:length(shared_cells{lvl, c}))});
            current_occupied_voxel_normals = cat(1, current_occupied_voxel_normals, occupied_voxel_normals{lvl, shared_cells{lvl, c}(2:length(shared_cells{lvl, c}))});
        end
        %Find the nearest neighbour to the current unique corner coordinate 
        diff = unique_coords{lvl}(c, :) - current_occupied_voxel_coords; 
        %euclid_dists = abs(arrayfun(@(idx) norm(diff(idx, :)), 1:size(diff, 1)));
        %Don't need "abs" below, because the diff values are squared
        euclid_dists = sqrt(sum(diff.^2, 2));   %This operation is faster than "arrayfun", above
        [min_euclid_dist{lvl}(c), nearest_voxel_ind] = min(euclid_dists);
        %Check if more than one same minimum Euclidean distance was found
        %above: if yes, this means that the current corner has more than
        %one equidistant nearest voxel. Currently do this only for the leaf
        %(voxel) level, where it is most likely to happen.
        if lvl == b + 1
            eqd_vox = find(euclid_dists == min(euclid_dists));
            if length(eqd_vox) > 1
                equidistant_voxel_inds{c} = eqd_vox;
            end
        end
        
        %Below is for debugging only
%         test_dists = find(euclid_dists == min_euclid_dist{lvl}(c));
%         if length(test_dists) > 1
%             test_dotproducts = zeros(length(test_dists), length(test_dists));
%             test_signs = zeros(length(test_dists), length(test_dists));
%             test_signs_final = zeros(length(test_dists), 1);
%             for s1 = 1:length(test_dists)
%                 for s2 = 1:length(test_dists)
%                     test_dotproducts(s1, s2) = dot(current_occupied_voxel_normals(test_dists(s1), :), current_occupied_voxel_normals(test_dists(s2), :));
%                     test_signs(s1, s2) = sign(test_dotproducts(s1, s2));
%                 end
%                 %Check if all the dot products are positive or not
%                 if (sum(test_signs(s1, :)) == length(test_dists))
%                     test_signs_final(s1) = 1;
%                 else
%                     test_signs_final(s1) = 0;
%                 end
%             end
%             if sum(test_signs_final) == length(test_dists)
%                 disp(['Found >1 equidistant nearest point for corner ' num2str(c) ' (normal sign difference: NO)']);
%             else
%                 disp(['Found >1 equidistant nearest point for corner ' num2str(c) ' (normal sign difference: YES)']);
%             end
%             disp('Normals:');
%             disp(num2str(current_occupied_voxel_normals(test_dists(:), :)));
%             disp('Dot products between normals:');
%             disp(num2str(test_dotproducts));
%             disp('Signs of dot products:');
%             disp(num2str(test_signs));
%             disp('Final signs (1 = row of signs sums to total no. of rows; 0 = row of signs sums to < total no. of rows)');
%             disp(num2str(test_signs_final));
%             disp(['Sum of final signs: ' num2str(sum(test_signs_final))]);
%         end

        %if (lvl < b + 1)
            %Store the (x, y, z) coordinates of the nearest voxel, for future
            %reference
            nearest_voxels{lvl}(c, 1:3) = current_occupied_voxel_coords(nearest_voxel_ind, :);
            %Get the normal for the nearest voxel (the normal is an (x, y, z)
            %triplet)
            normal_nearest_vox{lvl}(c, 1:3) = current_occupied_voxel_normals(nearest_voxel_ind, :); 
            %For the same corner point at the prevoius, lower octree level, 
            %check if a smaller min_euclid_dist (absolute value of the control 
            %point) is found than the current min_euclid_dist; keep the 
            %smallest min_euclid_dist for this corner. Update the nearest_voxel
            %coordinates and normal for this control point, accordingly.
            if lvl > 1
                %Initialize the true_min (true minimum control point value) to
                %the current absolute control point value
                %true_min = min_euclid_dist;
                true_min = min_euclid_dist{lvl}(c);
                %Currently, we are only considering one lower octree level
                lower_lvls = lvl - 1;            
                %Find the first row index where the current corner coordinates
                %are the same as the (x, y, z) coordinates of a corner at 
                %lower_lvls. We can stop after finding the first row index, 
                %because any other row indices that could be found would 
                %correspond to the same corner coordinate and therefore the 
                %same control point.
                corner_index = find(sum(abs(unique_coords{lvl}(c, :) - corner_coords{lower_lvls}), 2) == 0, 1);    
                ctrlpt_index = ctrl_pts_pointers{lower_lvls}(corner_index);
                if (~isempty(corner_index))&&(~isempty(control_points{lower_lvls}))
                    %All the indices in corner_index should be the same, so
                    %just pick the first one
                    %ctrlpt_index = ctrl_pts_pointers{lower_lvls}(corner_index(1));
                    %Find the control point value corresponding to the above 
                    %ctrlpt_index
                    test_ctrlpt = control_points{lower_lvls}(ctrlpt_index);
                    %Check if the absolute control point value found above is
                    %smaller than the current true_min: if yes, change the
                    %value of true_min to this new value and update the nearest
                    %voxel coordinates and normal accordingly.
                    if abs(test_ctrlpt) < true_min
                        true_min = abs(test_ctrlpt);
                        nearest_voxels{lvl}(c, 1:3) = nearest_voxels{lower_lvls}(ctrlpt_index, 1:3);
                        normal_nearest_vox{lvl}(c, 1:3) = normal_nearest_vox{lower_lvls}(ctrlpt_index, 1:3);
                    end
                end
                min_euclid_dist{lvl}(c) = true_min;
            end %End check if lvl > 1
        %If we are at the voxel level
        %elseif (lvl == b + 1)
%             %If, for the current corner, only one nearest voxel (midpoint) 
%             %has been found
%             if isempty(equidistant_voxel_inds{c})
%                 %Get the (x, y, z) coordinates of the nearest voxel
%                 %(midpoint)
%                 nearest_voxels{lvl}(c, 1:3) = current_occupied_voxel_coords(nearest_voxel_ind, :);
%                 %Get the normal for the nearest voxel (the normal is also
%                 %an (x, y, z) triplet)
%                 normal_nearest_vox{lvl}(c, 1:3) = current_occupied_voxel_normals(nearest_voxel_ind, :);       
%                 %Compute the difference vector between the current corner's
%                 %(x, y, z) coordinates and the nearest voxel's (x, y, z)
%                 %coordinates
%                 difference_vectors{lvl}(c, 1:3) = unique_coords{lvl}(c, :) - nearest_voxels{lvl}(c, 1:3);
%                 %Compute the dot product between the difference vector and 
%                 %normal vector for the current corner, found above
%                 dot_products{lvl}(c) = dot(difference_vectors{lvl}(c, 1:3), normal_nearest_vox{lvl}(c, 1:3), 2); 
%                 %To obtain the signed distance value for the control point at
%                 %the current corner, take the sign of the dot product computed
%                 %above for the current corner, and put it in front of the 
%                 %corresponding minimum Euclidean distance for this corner. A 
%                 %negative dot product indicates that the nearest voxel normal 
%                 %is pointing AWAY from the corresponding corner, so this corner 
%                 %is INSIDE the surface implied by the point cloud in that 
%                 %voxel. A positive dot product indicates that the nearest voxel 
%                 %normal is pointing TOWARDS the corresponding corner, so this 
%                 %corner is OUTSIDE the surface of the point cloud. 
%                 control_points{lvl}(c) = sign(dot_products{lvl}(c)).*(min(euclid_dists))';
%             %If, for the current corner, more than one equidistant voxel
%             %(midpoint) has been found
%             else
%                 %Get the (x, y, z) coordinates of all the equidistant 
%                 %nearest voxels (i.e., their midpoints)
%                 nearest_vox_coords = current_occupied_voxel_coords(equidistant_voxel_inds{c}, :);
%                 %Get the normal vector for each of the equidistant nearest
%                 %voxels
%                 normals_nearest_vox = current_occupied_voxel_normals(equidistant_voxel_inds{c}, :);
%                 %Compute the difference vector between the current corner's
%                 %(x, y, z) coordinates and each of the voxel coordinates in
%                 %nearest_vox_coords
%                 diff_vecs = unique_coords{lvl}(c, :) - nearest_vox_coords;
%                 %Compute the dot product between each difference vector and 
%                 %the corresponding nearest voxel's normal vector
%                 dot_prods = dot(diff_vecs, normals_nearest_vox, 2);
%                 %Find the average dot product from dot_prods
%                 mean_dp = mean(dot_prods);
%                 dot_products{lvl}(c) = mean_dp;
%                 %Set the control point for the current corner as the
%                 %average of the dot products found above
%                 control_points{lvl}(c) = dot_products{lvl}(c);               
% %                 %Compute the average of each component (x, y, z) of
% %                 %nearest_vox_coords
% %                 nearest_voxels{lvl}(c, 1:3) = mean(nearest_vox_coords);
% %                 %Compute the difference vector between the current corner's
% %                 %(x, y, z) coordinates and the average nearest voxel
% %                 %coordinates
% %                 difference_vectors{lvl}(c, 1:3) = unique_coords{lvl}(c, :) - nearest_voxels{lvl}(c, 1:3);
%   
% %                 %Average the normals in normals_nearest_vox
% %                 normal_nearest_vox{lvl}(c, 1:3) = mean(normals_nearest_vox);
% %                 %Compute the difference vector between the current corner's
% %                 %(x, y, z) coordinates and each of the equdistant nearest
% %                 %voxel's (x, y, z) coordinates
% %                 diff_vecs = unique_coords{lvl}(c, :) - nearest_vox_coords;
% %                 %Average the signs of each of the (x, y, z) components of 
% %                 %the difference vectors computed above, but leave the
% %                 %absolute values of the components the same (so, each 
% %                 %difference vector should have components of +/- 0.5)
% %                 difference_vectors{lvl}(c, 1:3) = 0.5.*sign(mean(sign(diff_vecs)));
%             end %End check if(isempty(equidistant_voxel_inds{c}))

        %end %End check if (lvl < b + 1)
        disp(['Finished corner ' num2str(c) '/' num2str(size(unique_coords{lvl}, 1))]);
    end %End c (current unique coordinate)
    %dot_products{lvl} = dot_products{lvl}'; %Make column vector
    %control_points{lvl} = control_points{lvl}'; %Make column vector
    %if lvl < b + 1
        %Compute the difference vector between each unique corner point's 
        %(x, y, z) coordinates and the corresponding nearest voxel's (x, y, z) 
        %coordinates, at the current octree level
        difference_vectors{lvl} = unique_coords{lvl} - nearest_voxels{lvl};
        %Compute the dot product between each difference vector and the
        %corresponding nearest voxel's normal vector, at the current octree
        %level
        dot_products{lvl} = dot(difference_vectors{lvl}, normal_nearest_vox{lvl}, 2);
        %To obtain the signed distance value for each control point at the
        %current octree level, take the sign of each dot product computed above
        %and add it to the corresponding min_euclid_dist. A negative dot
        %product indicates that the nearest voxel normal is pointing AWAY from
        %the corresponding corner, so this corner is INSIDE the surface
        %implied by the point cloud in that octree cell. A positive dot product
        %indicates that the nearest voxel normal is pointing TOWARDS the
        %corresponding corner, so this corner is OUTSIDE the surface of the 
        %point cloud. 
        control_points{lvl} = sign(dot_products{lvl}).*min_euclid_dist{lvl}';   
        disp('------------------------------------------------------------');
    %end
end %End lvl
ctrlpts_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all control points: ' num2str(ctrlpts_time) ' seconds']);
disp('************************************************************');

%For debugging purposes, print out how many voxels, if any, end up having 
%all 8 of their control points with the same sign (+/-) ...
[~, ptcloud, ~] = plyRead(ptcloud_file);
%Accumulate all the control points (not just the unique ones) into one cell
%array
all_ctrlpts = get_all_ctrlpts(control_points, ctrl_pts_pointers); 
same_sign_cntr = 0;
vox = 0;
same_sign_voxels = [];
disp(' ');
for i = 1:8:(length(all_ctrlpts{max_lvl}) - 7) 
    vox = vox + 1;
    %Get all 8 control points for the corners of the current voxel
    current_ctrlpts = all_ctrlpts{max_lvl}(i:(i + 7));
    %Check the signs of current_ctrlpts (assume here that none of the
    %control points have a value of 0)
    if (sum(sign(current_ctrlpts)) == length(current_ctrlpts))||(sum(sign(current_ctrlpts)) == -length(current_ctrlpts))
        same_sign_cntr = same_sign_cntr + 1;
        disp(' ');
        disp(['Voxel ' num2str(vox) ' (' num2str(ptcloud(vox, 1:3)) ') has all control points with the same sign']);
        disp('Voxel normal:');
        disp(normals_sorted(vox, :));
        %Store this voxel, for debugging purposes
        same_sign_voxels((size(same_sign_voxels, 1) + 1), 1:3) = ptcloud(vox, 1:3);
        %Get the corner coordinates for each corner of this voxel
        vox_corners = corner_coords{max_lvl}((i:(i + 7)), 1:3);
        disp('Voxel corner coordinates:');
        disp(num2str(vox_corners));
        %Get the nearest voxel midpoint found for each corner of the 
        %current voxel
        curr_nearest_vox = nearest_voxels{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + 7)), 1:3);
        disp('Nearest voxel midpoint found for each corner of the current voxel: ');
        disp(num2str(curr_nearest_vox));
        %Get the difference vector for each of the 8 corners of this voxel
        %(difference vector between each corner and the chosen voxel 
        %midpoint)
        vox_diff_vecs = difference_vectors{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + 7)), 1:3);
        disp('Voxel corner difference vectors:');
        disp(num2str(vox_diff_vecs));
        %Get the dot product between the difference vector of each corner
        %of this voxel and the normal vector of the chosen nearest voxel
        %midpoint (this may not be the midpoint of the current voxel, as
        %neighbouring voxels' midpoints are the same distance away)
        vox_dp = dot_products{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + 7)));
        disp('Voxel corner dot products:');
        disp(num2str(vox_dp));
    end
end
disp(' ');
disp(['TOTAL number of voxels with all control points having the same sign: ' num2str(same_sign_cntr) '/' num2str(length(all_ctrlpts{max_lvl})/8) ' (' num2str((same_sign_cntr/(length(all_ctrlpts{max_lvl})/8))*100) '%)']);
%Plot voxels that have the same control point signs
if same_sign_cntr > 0
    figure;
    scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, 'filled', 'MarkerFaceColor', 'b');
    hold on;
    scatter3(same_sign_voxels(:, 1), same_sign_voxels(:, 2), same_sign_voxels(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
    axis equal; axis off;
    title('Voxels with Same-Sign Control Points at Encoder');
    legend('Voxels with Different-Sign Control Points', 'Voxels with Same-Sign Control Points', 'Location', 'best');
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

%Pick a level from which to start the wavelet analysis (can go from the
%root (start_lvl = 1) up to 2 levels before the leaves (start_lvl = b - 1))
%start_lvl = 1;
%Choose a quantization stepsize, for uniform scalar quantization of the
%control points at the chosen base level (start_lvl) and all the wavelet
%coefficients computed at every other level
%q_stepsize = 1;
%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) across all 
%octree blocks and levels, starting from start_lvl and going up to one
%level before the leaves
%wavelet_coeffs = cell((b + 1 - start_lvl), 1);  %These will be quantized coefficients
%wavelet_coeffs = cell(b, 1);  %These will be quantized coefficients
wavelet_coeffs = cell((b + 1), 1);  %These will be quantized coefficients
%Initialize a cell array to store the reconstructed signal at each vertex
%of each octree cell at every level from start_lvl to one level before the
%leaves
reconstructed_control_points = cell(size(control_points, 1), 1);
%Initialize a flag that indicates whether or not we want to display
%visualizations at a given octree level (0 => no; 1 => yes). IMPORTANT:
%This must always be initialized to 0; it will be automatically set to 1
%later, if the vis_levels_wavelet array (see below) is not empty.
w_vis_flag = 0;
%Octree level(s) for which we want to visualize the wavelet analysis steps
%(leave the below array empty if you do not wish to visualize any wavelet 
%analysis steps) 
vis_levels_wavelet = [];    

if ~isempty(vis_levels_wavelet)
    %Read in the input point cloud (assume PLY format)
    [~, ptcloud, ~] = plyRead(ptcloud_file);
end

%Quantize all the control points at octree level start_lvl
for cpt = 1:size(control_points{start_lvl}, 1)
    control_points{start_lvl}(cpt, :) = quantize_uniform_scalar(control_points{start_lvl}(cpt, :), q_stepsize);
end
%Add these control points to reconstructed_control_points
reconstructed_control_points{start_lvl} = control_points{start_lvl};

%For each octree level, starting from start_lvl and working up to 1 level
%before the leaves ...
tic;
%for lvl = start_lvl:b
for lvl = start_lvl:(max_lvl - 1)
    disp(['Computing wavelet coefficients between octree levels ' num2str(lvl) ' and ' num2str(lvl + 1) ' ...']);
    disp('------------------------------------------------------------');
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Extract the cell's 8 corner coordinates. This cell will represent
        %our parent cell at the current level, since we will compute
        %wavelet coefficients for its child cells.
        parent_corner_coords = corner_coords{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);
        %Get the pointer to the first child cell of the current parent cell
        child_ptr = myOT.FirstChildPtr{lvl}(occ_cell);
        %For each child cell of the current cell ...
        for child = 1:myOT.ChildCount{lvl}(occ_cell)
            %Get the 8 corner coordinates of the current child cell. NOTE:
            %Below, "child" is converted to type double, because it is
            %uint8 by default, and if the result of any of the 
            %multiplications is > 255, the answer will be truncated and the
            %final result will be incorrect.
            child_corner_coords = corner_coords{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
            %For each corner of the current child cell ...
            for cnr = 1:8
                %Check if the current corner already has a wavelet 
                %coefficient associated with it (since some corners will be 
                %shared amongst different octree cells and we want to make
                %sure that we process only the UNIQUE corner vertices at 
                %each octree level) 
                possible_inds = ctrl_pts_pointers{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
                if ~isempty(wavelet_coeffs{lvl + 1})
                    if length(wavelet_coeffs{lvl + 1}) >= possible_inds(cnr) 
                        if ~isempty(wavelet_coeffs{lvl + 1}(possible_inds(cnr)))
                        %if (wavelet_coeffs{lvl + 1}(possible_inds(cnr))) ~= 0
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
                %Display visualizations only for the chosen octree level(s)
                if ~isempty(find(lvl == vis_levels_wavelet))
                    w_vis_flag = 1;
                    figure;
                    %Plot the input point cloud
                    scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
                    axis equal; axis off;
                    hold on;
                    %Outline the current parent cell in thick black lines 
                    %(i.e., connect all the vertices in parent_corner_coords  
                    %to the other vertices in parent_corner_coords, with 
                    %which they share an edge). Since we are only plotting
                    %one cell, we know in advance how the vertices are 
                    %connected, so first construct a matrix of edge indices.
                    edges_onecell = zeros(8, 3);    %Each corner vertex is connected to 3 others
                    edges_onecell(1, :) = [2, 4, 5];    %Corner 1 is connected to corners 2, 4, and 5
                    edges_onecell(2, :) = [1, 3, 6];
                    edges_onecell(3, :) = [2, 4, 7];
                    edges_onecell(4, :) = [1, 3, 8];
                    edges_onecell(5, :) = [1, 6, 8];
                    edges_onecell(6, :) = [2, 5, 7];
                    edges_onecell(7, :) = [3, 6, 8];
                    edges_onecell(8, :) = [4, 5, 7];
                    %Plot the root octree cell, for reference
                    for rootv1 = 1:8
                        for rootv2 = 1:3
                            h_w(1) = plot3([corner_coords{1}(rootv1, 1) corner_coords{1}(edges_onecell(rootv1, rootv2), 1)], [corner_coords{1}(rootv1, 2) corner_coords{1}(edges_onecell(rootv1, rootv2), 2)], [corner_coords{1}(rootv1, 3) corner_coords{1}(edges_onecell(rootv1, rootv2), 3)], 'r');
                            hold on;
                        end
                    end
                    %Connect the vertices in the current parent cell, 
                    %according to edges_onecell
                    for pvtx1 = 1:8
                        for pvtx2 = 1:3
                            h_w(2) = plot3([parent_corner_coords(pvtx1, 1) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 1)], [parent_corner_coords(pvtx1, 2) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 2)], [parent_corner_coords(pvtx1, 3) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 3)], 'k', 'LineWidth', 3);
                            hold on;
                        end
                    end
                    %Outline the current child cell in thick, dotted black 
                    %lines, by connecting the vertices in child_corner_coords 
                    %according to edges_onecell
                    for cvtx1 = 1:8
                        for cvtx2 = 1:3
                            h_w(3) = plot3([child_corner_coords(cvtx1, 1) child_corner_coords(edges_onecell(cvtx1, cvtx2), 1)], [child_corner_coords(cvtx1, 2) child_corner_coords(edges_onecell(cvtx1, cvtx2), 2)], [child_corner_coords(cvtx1, 3) child_corner_coords(edges_onecell(cvtx1, cvtx2), 3)], '--k', 'LineWidth', 3);
                            hold on;
                        end
                    end               
                    %Plot the current corner point (in red, with a black 
                    %outline) on the current figure
                    h_w(4) = scatter3(child_corner_coords(cnr, 1), child_corner_coords(cnr, 2), child_corner_coords(cnr, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
                    hold on;   
                end %End first part of the visualization code for the current corner     
                
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
                    %Get the row indices of the parent vertices on this 
                    %face
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
                    %parent_cnr_index = find(ismember(parent_corner_coords, child_corner_coords(cnr, :), 'rows') > 0);
                    %Find the row index of this coordinate inside 
                    %corner_coords at the parent octree level
                    %parent_row_index = find(ismember(corner_coords{lvl}, parent_corner_coords(parent_cnr_index, :), 'rows') > 0);
                    %parent_row_index = find(ismember(corner_coords{lvl}, child_corner_coords(cnr, :), 'rows') > 0); 
                    parent_row_index = find(sum(abs(child_corner_coords(cnr, :) - corner_coords{lvl}), 2) == 0, 1); %Faster than "ismember"
                    
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
                    %except insert a 0 here in the wavelet_coeffs cell 
                    %array, and transfer the reconstructed control point 
                    %over from the parent corner at the previous octree 
                    %level (we want the reconstructed_control_points values  
                    %and wavelet_coeffs values at each level to correspond
                    %to the same locations in unique_coords)
                    wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = 0;
                    reconstructed_control_points{lvl + 1}(possible_inds(cnr)) = reconstructed_control_points{lvl}(parent_ctrlpt_index);
                    if w_vis_flag == 1
                        hold off;
                        legend(h_w, 'Root cell (for reference)', ['Current parent cell (at level ' num2str(lvl) ')'], ['Current child cell (at level ' num2str(lvl + 1) ')'], 'Current child corner', 'Location', 'best');
                        title({'Computing Wavelet Coefficients', ['for Each Child of Each Occupied Octree Cell at Level ' num2str(lvl)]});
                        w_vis_flag = 0; %Reset flag for next octree cell
                    end
                    continue;
                %If this corner vertex lies on a parent edge
                elseif on_p_edge == 1
                    %Get the indices of the rows of these coordinates in
                    %parent_corner_coords (should be only 2 rows) 
                    if w_vis_flag == 1
                        for pc = 1:length(parent_row_inds)
                            %Circle the parent corners (in blue)
                            h_w(5) = scatter3(parent_corner_coords(parent_row_inds(pc), 1), parent_corner_coords(parent_row_inds(pc), 2), parent_corner_coords(parent_row_inds(pc), 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                            hold on;
                        end
                    end
                    %Get the Bezier control points stored at the corner 
                    %vertices defined by parent_row_inds
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt1 = reconstructed_control_points{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstructed_control_points{lvl}(ctrlpt2_ptr);
                    %Average the signal (Bezier control points) on the 2
                    %vertices of the parent edge
                    avg_signal = (ctrlpt1 + ctrlpt2)/2;
                %If this corner vertex lies on a parent face 
                elseif on_p_face == 1
                    %Get the indices of the rows of these coordinates in
                    %parent_corner_coords (should be 4 rows) 
                    %parent_row_inds = find(sum(ismember(parent_corner_coords, cell_corner_coords(cnr, :)), 2) == 1);
                    if w_vis_flag == 1
                        for pc = 1:length(parent_row_inds)
                            %Circle the parent corners (in blue)
                            h_w(5) = scatter3(parent_corner_coords(parent_row_inds(pc), 1), parent_corner_coords(parent_row_inds(pc), 2), parent_corner_coords(parent_row_inds(pc), 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                            hold on;
                        end
                    end
                    %Get the Bezier control points stored at the corner
                    %vertices defined by parent_row_inds
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt3_ptr = all_ctrlpts_ptrs(parent_row_inds(3));
                    ctrlpt4_ptr = all_ctrlpts_ptrs(parent_row_inds(4));
                    ctrlpt1 = reconstructed_control_points{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstructed_control_points{lvl}(ctrlpt2_ptr);
                    ctrlpt3 = reconstructed_control_points{lvl}(ctrlpt3_ptr);
                    ctrlpt4 = reconstructed_control_points{lvl}(ctrlpt4_ptr);
                    %Average the signal (Bezier control points) on the 4 
                    %vertices of the parent face
                    avg_signal = (ctrlpt1 + ctrlpt2 + ctrlpt3 + ctrlpt4)/4;  
                %If this corner vertex lies somewhere in the centre of the
                %parent's block (i.e., neither on a parent's edge or on a
                %parent's face)
                else                
                    if w_vis_flag == 1
                        for pc = 1:8
                            %Circle the parent corners (in blue)
                            h_w(5) = scatter3(parent_corner_coords(pc, 1), parent_corner_coords(pc, 2), parent_corner_coords(pc, 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                            hold on;
                        end
                    end
                    %Get the Bezier control points stored at each of the 8
                    %corner vertices of the parent cell
                    ctrlpts_pointers = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpts = reconstructed_control_points{lvl}(ctrlpts_pointers);
                    %Average the signal (Bezier control points) on the 8
                    %vertices of the parent block
                    avg_signal = mean(ctrlpts);
                end
                %Get the Bezier control point stored at the current corner
                %of the current child cell
                child_ctrlpt = control_points{lvl + 1}(possible_inds(cnr));
                %Subtract the average signal from the signal (Bezier
                %control point) at the current child vertex. The result is 
                %the high-pass transform (wavelet) coefficient of the child 
                %vertex.
                wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = child_ctrlpt - avg_signal;
                %Quantize the wavelet coefficient computed above
                wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(possible_inds(cnr)), q_stepsize);
                %Add the quantized wavelet coefficient to avg_signal, to 
                %obtain the reconstructed signal (control point) at the 
                %current child vertex 
                reconstructed_control_points{lvl + 1}(possible_inds(cnr)) = wavelet_coeffs{lvl + 1}(possible_inds(cnr)) + avg_signal;
                if w_vis_flag == 1
                    hold off;
                    legend(h_w, 'Root cell (for reference)', ['Current parent cell (at level ' num2str(lvl) ')'], ['Current child cell (at level ' num2str(lvl + 1) ')'], 'Current child corner', 'Parent corners', 'Location', 'best');
                    title({'Computing Wavelet Coefficients', ['for Each Child of Each Occupied Octree Cell at Level ' num2str(lvl)]});
                    w_vis_flag = 0; %Reset flag for next octree cell
                end   
            end %End corners
        end %End children
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
    end %End occupied cells at current level
    %Arrange all the wavelet coefficients produced for the child octree 
    %level, into a column vector instead of a row vector (purely for 
    %visualization reasons: it is easier to scroll through a long column 
    %vector than a row vector). Also, wavelet_coeffs will then be in the 
    %same format as the control_points cell array.
    wavelet_coeffs{lvl + 1} = (wavelet_coeffs{lvl + 1})';
    %Also arrange all reconstructed control points at the next level, into
    %a column vector
    reconstructed_control_points{lvl + 1} = (reconstructed_control_points{lvl + 1})';
end %End current octree level
waveletcfs_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all wavelet coefficients: ' num2str(waveletcfs_time) ' seconds']);
disp('************************************************************');

%Plot the distribution of reconstructed control points (low-pass
%coefficietns) at different octree levels
for lvl = start_lvl:max_lvl
    figure;
    histogram(reconstructed_control_points{lvl});
    title(['Histogram of (Encoder-Reconstructed) Control Points at Octree Level ' num2str(lvl)]);
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_ctrlpts_lvl' num2str(lvl)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_ctrlpts_lvl' num2str(lvl)], '-dpdf');
end
%Plot the distribution of quantized wavelet coefficients at different
%octree levels
for lvl = start_lvl:(max_lvl - 1)
    figure;
    histogram(wavelet_coeffs{lvl + 1});
    title(['Histogram of Quantized Wavelet Coefficients between Octree Levels ' num2str(lvl) ' and ' num2str(lvl + 1)]);
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_wavcfs_lvls' num2str(lvl) '-' num2str(lvl + 1)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\histogram_wavcfs_lvls' num2str(lvl) '-' num2str(lvl + 1)], '-dpdf');
end

%---------------- Checking for Zero Wavelet Coefficients -----------------%

disp(' ');
disp('------- Checking for All Zero Wavelet Coefficients  --------');
disp(' ');

%Cell array to store the indices of the occupied octree cells at each level
%that have zero wavelet coefficients on all of their corners
all_zero_wav_cfs = cell(size(wavelet_coeffs));
%Cell array to store the midpoints of the occupied octree cells at each
%level that have zero wavelet coefficients on all of their corners (only
%needed for display purposes, later)
all_zero_cell_midpoints = cell(size(wavelet_coeffs));
%Cell array to store the indices of ctrl_pts_pointers for the occupied
%octree cells at each level that have zero wavelet coefficients on all of
%their corners (only needed for display purposes, later)
all_zero_cell_ctrlpts_ptrs = cell(size(wavelet_coeffs));

%At each octree level, check which occupied octree cells (if any) have zero
%wavelet coefficients on all of their corners
for lvl = (start_lvl + 1):max_lvl
    disp(['Processing octree level ' num2str(lvl) ' ...']);
    %Counter for number of occupied cells at this level, which contain all
    %zero wavelet coefficients
    zw_cntr = 1;
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Get the wavelet coefficient for each of the 8 corners of this cell
        current_wavelet_coeffs = wavelet_coeffs{lvl}(ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8)));
        %If the wavelet coefficients at all the corners of this cell are 0
        if isempty(find(current_wavelet_coeffs ~= 0))
            disp(['All zero wavelet coefficients for occupied cell ' num2str(occ_cell)]);
            all_zero_wav_cfs{lvl}(zw_cntr) = occ_cell; 
            %Get the midpoints of the current cell (only needed for display
            %purposes, later)
            all_zero_cell_midpoints{lvl}(zw_cntr, 1:3) = mean(corner_coords{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :), 1);
            %Store the control points pointers for each corner of this cell
            %(only needed for display purposes, later)
            all_zero_cell_ctrlpts_ptrs{lvl}(zw_cntr, 1:8) = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8));
            zw_cntr = zw_cntr + 1;
        end
    end
    disp(['TOTAL no. of occupied cells with all zero wavelet coefficients at this level: ' num2str(length(all_zero_wav_cfs{lvl}))]);
    disp('------------------------------------------------------------');
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

%Determine if octree and wavelet coefficient tree pruning is needed
if prune_flag == 1
    %-------------------------- Pruning Octree ---------------------------%

    disp(' ');
    disp('-------------------- Pruning Octree ------------------------');
    disp(' ');

    %Create a copy of the myOT.OccupancyCode cell array, which will contain
    %the pruned set of occupancy codes
    pruned_occupancy_codes = myOT.OccupancyCode;
    %Create a cell array that will contain the indices of the octree cells 
    %whose occupancy codes (which describe the locations of their occupied 
    %children) are to be pruned away at each octree level, from
    %pruned_occupancy_codes. These will be the internal cells only, not the
    %leaves.
    toprune = cell(size(myOT.OccupancyCode));
    %Create a supplementary cell array that will contain a 1 where the 
    %corresponding occupied cell is a leaf node (after pruning) and a 0 
    %where the corresponding occupied cell is internal (not a leaf node).
    post_pruning_array = cell(size(myOT.OccupancyCode));
    %For each octree level that contains wavelet coefficients, except the 
    %voxel level (since this level does not have its own occupancy codes,
    %as the voxels do not have children), initialize all locations in 
    %post_pruning_array to 0 (the leaf cells ("1" locations) will be
    %determined after pruning)
    for lvl = (start_lvl + 1):(max_lvl - 1) %Here assume max_lvl = b + 1
        post_pruning_array{lvl} = zeros(length(myOT.OccupancyCode{lvl}), 1);
    end   
    %Create a cell array that will contain the indices of the octree blocks
    %that need to be pruned away from post_pruning_array, to account for 
    %the internal nodes that were removed in pruned_occupancy_codes.
    %NOTE: pruned_occupancy_codes and post_pruning_array will not be the 
    %same size at the end: pruned_occupancy_codes will contain one 
    %occupancy code per INTERNAL octree cell at each octree level, whereas 
    %post_pruning_array will contain a 1 or 0 bit for EACH of the remaining 
    %octree cells after pruning, including the internal cells AND the 
    %leaves (except the leaves at level b + 1).
    toprune2 = cell(size(myOT.OccupancyCode));

%     %For each octree level that contains wavelet coefficients, except the 
%     %leaves ...
%     for lvl = (start_lvl + 1):(max_lvl - 1) %Here assume max_lvl = b + 1
%         %If this level contains any occupied cells that have all zero
%         %wavelet coefficients ...
%         if ~isempty(all_zero_wav_cfs{lvl})
%             %In post_pruning_array at the current level, put a 1 in the
%             %location(s) corresponding to the occupied cell(s) found above,
%             %and a 0 in the other locations (1 indicates that the 
%             %corresponding cells are leaf cells after pruning, and 0 
%             %indicates that the corresponding cells are internal)
%             post_pruning_array{lvl} = zeros(length(myOT.OccupancyCode{lvl}), 1);
%             post_pruning_array{lvl}(all_zero_wav_cfs{lvl}) = 1;
%         end
%     end

    %The below code prunes descendants of octree cells that have all-zero
    %wavelet coefficients, without checking whether all these descendants
    %also have all-zero wavelet coefficients ...
    
%     %For each octree level that contains wavelet coefficients, except the
%     %leaves ...
%     for lvl = (start_lvl + 1):(max_lvl - 1) %Here assume max_lvl = b + 1
%         %If there exist any occupied octree cells at this level, which  
%         %contain all zero wavelet coefficients ...
%         if ~isempty(all_zero_wav_cfs{lvl})
%             %For each cell at this level, which contains all zero wavelet
%             %coefficients ...
%             for pruned_cell = all_zero_wav_cfs{lvl}
%                 %Add this cell's index to the toprune array
%                 toprune{lvl}(length(toprune{lvl}) + 1) = pruned_cell;
%                 if lvl < (max_lvl - 1)
%                     %Add the indices of each of this cell's children to 
%                     %toprune
%                     toprune{lvl + 1}((length(toprune{lvl + 1}) + 1):(myOT.ChildCount{lvl}(pruned_cell) + length(toprune{lvl + 1}))) = ((myOT.FirstChildPtr{lvl}(pruned_cell)):(myOT.FirstChildPtr{lvl}(pruned_cell) + uint32((myOT.ChildCount{lvl}(pruned_cell))) - 1));
%                     %Add the same indices to toprune2
%                     toprune2{lvl + 1}((length(toprune2{lvl + 1}) + 1):(myOT.ChildCount{lvl}(pruned_cell) + length(toprune2{lvl + 1}))) = ((myOT.FirstChildPtr{lvl}(pruned_cell)):(myOT.FirstChildPtr{lvl}(pruned_cell) + uint32((myOT.ChildCount{lvl}(pruned_cell))) - 1));
%                     %The current child cells will be the starting point for
%                     %looking for the other descendants of the current 
%                     %pruned_cell, at the next higher octree level
%                     parent_cells = ((myOT.FirstChildPtr{lvl}(pruned_cell)):(myOT.FirstChildPtr{lvl}(pruned_cell) + uint32((myOT.ChildCount{lvl}(pruned_cell))) - 1));
%                     %For each higher octree level, except the leaves (since
%                     %there are no occupancy codes at level b + 1 in 
%                     %myOT.OccupancyCode, as the leaf occupancies are 
%                     %described at level b)
%                     for lvl2 = (lvl + 1):(max_lvl - 2)
%                         %Get the descendants of the current pruned_cell, at 
%                         %level lvl2
%                         descendants = zeros(sum(myOT.ChildCount{lvl2}(parent_cells)), 1);
%                         first_children = myOT.FirstChildPtr{lvl2}(parent_cells);
%                         desc_cntr = 1;
%                         for fc = 1:length(first_children)
%                             descendants(desc_cntr:(desc_cntr + myOT.ChildCount{lvl2}(parent_cells(fc)) - 1)) = first_children(fc):(first_children(fc) + uint32(myOT.ChildCount{lvl2}(parent_cells(fc)) - 1));
%                             desc_cntr = desc_cntr + myOT.ChildCount{lvl2}(parent_cells(fc));
%                         end
%                         %Add the indices of descendants at lvl2 to toprune
%                         toprune{lvl2 + 1}((length(toprune{lvl2 + 1}) + 1):(length(descendants) + length(toprune{lvl2 + 1}))) = descendants;
%                         %Add the same indices to toprune2
%                         toprune2{lvl2 + 1}((length(toprune2{lvl2 + 1}) + 1):(length(descendants) + length(toprune2{lvl2 + 1}))) = descendants;
%                         %The new parent_cells, at the next octree level 
%                         %(lvl2 + 1), will be the current "descendants"
%                         parent_cells = descendants;
%                     end %End lvl2  
%                 end %End check if lvl < (max_lvl - 1)
%             end %End pruned_cell    
%             %Keep only the unique indices in toprune{lvl} and toprune2{lvl}
%             toprune{lvl} = unique(toprune{lvl}, 'stable');
%             toprune2{lvl} = unique(toprune2{lvl}, 'stable');
%             %Make toprune{lvl} and toprune2{lvl} column vectors rather than
%             %row vectors (only because they are easier to scroll through)
%             toprune{lvl} = toprune{lvl}';
%             toprune2{lvl} = toprune2{lvl}';
%             disp(['Level ' num2str(lvl) ' cells and their descendants have been marked for pruning']);
%             disp('------------------------------------------------------------');
%         end %End check if ~isempty(all_zero_wav_cfs{lvl})
%     end %End lvl

    %The below code prunes branches (occupancy codes) of octree cells that 
    %have all zero wavelet coefficients throughout the branch, so the 
    %pruning is done bottom-up from the voxel level ...
    
    %For each octree level that contains wavelet coefficients ...
    %for lvl = max_lvl:-1:(start_lvl + 2) %Here assume max_lvl = b + 1
    %Start from the voxel level ...
    for lvl = max_lvl %Here assume max_lvl = b + 1
        %If there exist any occupied octree cells at this level, which  
        %contain all zero wavelet coefficients ...
        if ~isempty(all_zero_wav_cfs{lvl})
            disp(['Found octree cells with all-zero wavelet coefficients at level ' num2str(lvl) ': ']);
            disp(num2str(all_zero_wav_cfs{lvl}'));
            %For each of these cells, check all of its ancestors to see if
            %there exist any branches of the octree that can be pruned off
            %(i.e., in the case that all of the octree cells in a branch
            %have all-zero wavelet coefficients)
            for az_cell = all_zero_wav_cfs{lvl}
                %Find the parent of the current az_cell
                parent_cell = myOT.ParentPtr{lvl}(az_cell);
                %Find all the children of parent_cell (one of these
                %children will be the current az_cell)
                children = ((myOT.FirstChildPtr{lvl - 1}(parent_cell)):(myOT.FirstChildPtr{lvl - 1}(parent_cell) + uint32((myOT.ChildCount{lvl - 1}(parent_cell))) - 1));
                %Check if this parent's other children also have all-zero
                %wavelet coefficients
                az_children = children(find(ismember(children, all_zero_wav_cfs{lvl}) == 1));
                %Check if parent_cell has all-zero wavelet coefficients
                if (ismember(parent_cell, all_zero_wav_cfs{lvl - 1}) == 1)
                    %If all of its children also have all-zero wavelet
                    %coefficients
                    if length(az_children) == length(children)
                        %Mark the occupancy code for the current 
                        %parent_cell for pruning (since its children will
                        %be pruned off)
                        toprune{lvl - 1}(length(toprune{lvl - 1}) + 1) = parent_cell;
                        %If children are not at the voxel level, mark the 
                        %occupancy codes of parent_cell's occupied children
                        %for pruning too (assuming here that the occupied 
                        %children will always have occupied children of 
                        %their own)
                        if (lvl < (b + 1))&&(~isempty(myOT.OccupancyCode{lvl}(children)))
                            toprune{lvl}((length(toprune{lvl}) + 1):(length(toprune{lvl}) + length(children))) = children;
                            %These children will not be leaves (since they
                            %will be pruned off), so we can mark them for 
                            %removal from post_pruning_array as well
                            toprune2{lvl}((length(toprune2{lvl}) + 1):(length(toprune2{lvl}) + length(children))) = children;
                        end
                        %Successively check each ancestor of the current
                        %parent_cell, to see if we can prune further up the 
                        %octree 
                        for lvl2 = (lvl - 1):-1:(start_lvl + 1)
                            parent2 = myOT.ParentPtr{lvl2}(parent_cell);
                            children2 = ((myOT.FirstChildPtr{lvl2 - 1}(parent2)):(myOT.FirstChildPtr{lvl2 - 1}(parent2) + uint32((myOT.ChildCount{lvl2 - 1}(parent2))) - 1));
                            az_children2 = children2(find(ismember(children2, all_zero_wav_cfs{lvl2}) == 1));
                            %If the ancestor parent has all-zero wavelet
                            %coefficients 
                            if (ismember(parent2, all_zero_wav_cfs{lvl2 - 1}) == 1)
                                %If its other children also have all-zero 
                                %wavelet coefficients
                                if(length(az_children2) == length(children2))
                                    %Mark the occupancy code for the
                                    %current ancestor parent for pruning
                                    toprune{lvl2 - 1}(length(toprune{lvl2 - 1}) + 1) = parent2;
                                    %Mark the current parent_cell for
                                    %pruning from post_pruning_array too, 
                                    %since it is not a leaf (as it will be
                                    %pruned off)
                                    toprune2{lvl2}(length(toprune2{lvl2}) + 1) = parent_cell;
                                    %If we are not at the voxel level, mark 
                                    %the occupancy codes of parent2's
                                    %occupied children for pruning too 
                                    %(assuming here that the occupied 
                                    %children will always have occupied 
                                    %children of their own)
                                    if (lvl2 < (b + 1))&&(~isempty(myOT.OccupancyCode{lvl2}(children2)))
                                        toprune{lvl2}((length(toprune{lvl2}) + 1):(length(toprune{lvl2}) + length(children2))) = children2;
                                        %These children will not be leaves 
                                        %(since they will be pruned off), 
                                        %so we can mark them for removal 
                                        %from post_pruning_array as well
                                        toprune2{lvl2}((length(toprune2{lvl2}) + 1):(length(toprune2{lvl2}) + length(children2))) = children2;
                                    end
                                    %Continue checking other ancestors,
                                    %with the current parent2 being the new 
                                    %parent_cell
                                    parent_cell = parent2;
                                    continue;  
                                %If the ancestor parent's children do not
                                %all have all-zero wavelet coefficients
                                else
                                    %We cannot prune the occupancy code of
                                    %the ancestor parent, but the current
                                    %parent_cell will become a leaf (its
                                    %occupancy code was already marked for
                                    %pruning, in the previous iteration of
                                    %lvl2)
                                    post_pruning_array{lvl2}(parent_cell) = 1;
                                    %Stop searching other ancestors; move
                                    %on to the next az_cell
                                    break;
                                end
                            end %End check if (ismember(parent2, all_zero_wav_cfs{lvl2 - 1}) == 1)
                            %If the ancestor parent does not have all-zero
                            %wavelet coefficients
                            if (ismember(parent2, all_zero_wav_cfs{lvl2 - 1}) == 0)
                                %But if all of its children do have
                                %all-zero wavelet coefficients
                                if (length(az_children2) == length(children2))
                                    %Mark the occupancy code of the current
                                    %ancestor parent, for pruning
                                    toprune{lvl2 - 1}(length(toprune{lvl2 - 1}) + 1) = parent2;
                                    %The ancestor parent will become a leaf
                                    post_pruning_array{lvl2 - 1}(parent2) = 1;
                                    %Mark the current parent_cell for
                                    %pruning from post_pruning_array, since
                                    %it is not a leaf
                                    toprune2{lvl2}(length(toprune2{lvl2}) + 1) = parent_cell;
                                    %If we are not at the voxel level, mark 
                                    %the occupancy codes of parent2's
                                    %occupied children for pruning too 
                                    %(assuming here that the occupied 
                                    %children will always have occupied 
                                    %children of their own)
                                    if (lvl2 < (b + 1))&&(~isempty(myOT.OccupancyCode{lvl2}(children2)))
                                        toprune{lvl2}((length(toprune{lvl2}) + 1):(length(toprune{lvl2}) + length(children2))) = children2;
                                        %These children will not be leaves
                                        %(since they will be pruned off),
                                        %so we can mark them for removal 
                                        %from post_pruning_array as well
                                        toprune2{lvl2}((length(toprune2{lvl2}) + 1):(length(toprune2{lvl2}) + length(children2))) = children2;
                                    end
                                    %Stop searching other ancestors; move
                                    %on to the next az_cell
                                    break;
                                %If the ancestor parent's children do not
                                %all have all-zero wavelet coefficients
                                else
                                    %The current parent_cell will become a
                                    %leaf
                                    post_pruning_array{lvl2}(parent_cell) = 1;
                                    %Stop searching other ancestors; move
                                    %on to the next az_cell
                                    break;
                                end
                            end %End check if (ismember(parent2, all_zero_wav_cfs{lvl2 - 1}) == 0)      
                        end %End lvl2
                    else
                        %Since all of the other children of parent_cell do
                        %not have all-zero wavelet coefficients, we cannot
                        %prune below parent_cell, so move on to checking
                        %the next az_cell
                        continue;
                    end %End check if parent_cell's other children have all-zero wavelet coefficients
                %If parent_cell does NOT have all-zero wavelet coefficients
                else
                    %But if all of its children do have all-zero wavelet
                    %coefficients
                    if length(az_children) == length(children)
                        %Mark the occupancy code of parent_cell for pruning
                        toprune{lvl - 1}(length(toprune{lvl - 1}) + 1) = parent_cell;
                        %The parent_cell will become a leaf after pruning
                        post_pruning_array{lvl - 1}(parent_cell) = 1;
                        %If we are not at the voxel level, mark the 
                        %occupancy codes of parent_cell's occupied children
                        %for pruning too (assuming here that the occupied 
                        %children will always have occupied children of 
                        %their own)
                        if (lvl < (b + 1))&&(~isempty(myOT.OccupancyCode{lvl}(children)))
                            toprune{lvl}((length(toprune{lvl}) + 1):(length(toprune{lvl}) + length(children))) = children;
                            %These children will not be leaves (since they
                            %will be pruned off), so we can mark them for 
                            %removal from post_pruning_array as well
                            toprune2{lvl}((length(toprune2{lvl}) + 1):(length(toprune2{lvl}) + length(children))) = children; 
                        end
                    else
                        %If neither the parent_cell nor all of its children
                        %have all-zero wavelet coefficients, we cannot
                        %prune here, so move on to checking the next
                        %az_cell
                        continue;
                    end                    
                end %End check if parent_cell has all-zero wavelet coefficients
            end %End az_cell
            
            if (lvl < (b + 1))&&(~isempty(toprune{lvl})) 
                %Keep only the unique indices in toprune{lvl}
                toprune{lvl} = unique(toprune{lvl}, 'stable');
                %Make toprune{lvl} a column vector rather than a row vector
                %(only because column vectors are easier to scroll through)
                toprune{lvl} = toprune{lvl}';
            end
            if (lvl < (b + 1))&&(~isempty(toprune2{lvl}))
                %Keep only the unique indices in toprune2{lvl}
                toprune2{lvl} = unique(toprune2{lvl}, 'stable');
                %Make toprune2{lvl} a column vector rather than a row 
                %vector (only because column vectors are easier to scroll
                %through)
                toprune2{lvl} = toprune2{lvl}';
            end
 
            disp(['Ancestors of Level ' num2str(lvl) ' have been checked for pruning']);
            disp('------------------------------------------------------------'); 
     
        end %End check if ~isempty(all_zero_wav_cfs{lvl})
    end %End lvl    
    
    %For each octree level covered in toprune ...
    for lvl = 1:size(toprune, 1)
        if ~isempty(toprune{lvl})
            %Prune the octree cells that have been marked
            pruned_occupancy_codes{lvl}(toprune{lvl}) = [];
            disp(['Level ' num2str(lvl) ' octree pruning done (removed occupancy codes of new leaf cells and their descendants)']);
            disp('------------------------------------------------------------');
        end
    end

    %For each octree level covered in toprune2 ...
    for lvl = 1:size(toprune2, 1)
        if ~isempty(toprune2{lvl})
            %Prune the post_pruning_array accordingly
            post_pruning_array{lvl}(toprune2{lvl}) = [];
            disp(['Level ' num2str(lvl) ' post-pruning array pruning done']);
            disp('------------------------------------------------------------');
        end
    end
    
    %----------------- Pruning Wavelet Coefficient Tree ------------------%

    disp(' ');
    disp('------------ Pruning Wavelet Coefficient Tree --------------');
    disp(' ');

    %Create a cell array that will contain indices of the wavelet 
    %coefficients that will be pruned off at each octree level
    wc_toprune = cell(size(wavelet_coeffs));
    %Create a copy of the wavelet_coeffs cell array, which will contain the 
    %pruned set of wavelet coefficients
    pruned_wavelet_coeffs = wavelet_coeffs;
    
    %For each octree level at which there are wavelet coefficients ...
    for lvl = (start_lvl + 1):max_lvl
        %Counter to keep track of the number of wavelet coefficients at the
        %current level, which will be pruned off
        wc_toprune_cntr = 1;
        %If we are not at the leaf level, we can look at the toprune2 cell
        %array (use this one instead of toprune, because toprune2 does not
        %mark leaf cells for pruning and we want to keep wavelet
        %coefficients for the leaf cells since their corner coordinates
        %will still need to be reconstructed at the decoder and we need to
        %have a wavelet coefficient per corner)
        if lvl < (b + 1)
            %If this level contains any octree cells whose occupancy 
            %codes have been pruned, except for leaf cells ...
            if ~isempty(toprune2{lvl})
                %Get the indices of these cells
                pruned_cells = toprune2{lvl};
                %For each unique corner at the current level ...
                for c = 1:length(wavelet_coeffs{lvl})
                    %Find the indices of the octree cells at this level, 
                    %which share this corner
                    current_shared_cells = shared_cells{lvl, c};
                    %Check if any of the shared cells found above are 
                    %inside pruned_cells: if ALL of the shared cells are 
                    %inside pruned_cells, then the wavelet coefficient for 
                    %this corner will be pruned off, so record its index
                    test = find(ismember(current_shared_cells, pruned_cells)) == 1;
                    if length(test) == length(current_shared_cells)
                        wc_toprune{lvl}(wc_toprune_cntr) = c;
                        wc_toprune_cntr = wc_toprune_cntr + 1;
                    end
                end %End c
            end %End check if ~isempty(toprune2{lvl})
        %If lvl == (b + 1)
        else  
            %The voxel level is not covered in the toprune2 array (since it
            %does not have occupancy codes of its own, as voxels have no
            %children), so for this level consider only the voxels that are 
            %descendants of octree cells that were pruned off at previous
            %octree levels. For each unique corner at the current level ...
            for c = 1:length(wavelet_coeffs{lvl})
                %Find the indices of the octree cells at this level,
                %which share this corner
                current_shared_cells = shared_cells{lvl, c};
                %Check if any of the current_shared_cells are children of
                %cells that were pruned off at the previous octree level:
                %if ALL of the shared cells are children of pruned-off
                %cells, then the wavelet coefficient for this corner will
                %be pruned off, so record its index
                parents = myOT.ParentPtr{lvl}(current_shared_cells);
                %Use toprune below, instead of toprune2, because we want to
                %consider leaf cells at the previous level as well as 
                %internal cells, but the leaf cells would not be marked in 
                %toprune2
                test = find(ismember(parents, toprune{lvl - 1})) == 1;
                if length(test) == length(parents)
                    wc_toprune{lvl}(wc_toprune_cntr) = c;
                    wc_toprune_cntr = wc_toprune_cntr + 1;
                end    
            end %End c    
        end %End check if lvl < (b + 1)
        %Make wc_toprune{lvl} a column vector rather than a row vector
        %(only because this is easier to scroll through)
        wc_toprune{lvl} = wc_toprune{lvl}';
        %Prune off all the wavelet coefficients whose indices are inside
        %wc_toprune at the current octree level
        if ~isempty(wc_toprune{lvl})
            pruned_wavelet_coeffs{lvl}(wc_toprune{lvl}) = [];
            disp(['Level ' num2str(lvl) ' done']);
            disp('------------------------------------------------------------');
        end
    end %End lvl   
end %End check if prune_flag == 1
 
%------------------------------- Encoding --------------------------------%

disp(' ');
disp('--------------- Encoding for Transmission ------------------');
disp(' ');

%Read in the input point cloud (assume PLY format), in order to get the
%total number of input voxels
[~, ptcloud, ~] = plyRead(ptcloud_file);

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
%Plot a histogram of the occupancy codes inside occ_codes_array
figure;
histogram(occ_codes_array);
title(['Histogram of Octree Occupancy Codes from Level 1-' num2str(max_lvl - 1)]);
%Save the above histogram as a MATLAB figure and as a PDF image in our
%network directory (NB: The '-bestfit' option maximizes the size of the 
%figure to fill the page, but preserves the aspect ratio of the figure. 
%The figure might not fill the entire page. This option leaves a 
%minimum page margin of .25 inches).
savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_occ_codes_lvl1-' num2str(max_lvl - 1)]);
print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_occ_codes_lvl1-' num2str(max_lvl - 1)], '-dpdf');
%Compute the number of bits required for these occupancy codes
bits_occ_codes_persymbol = entropy_calc(occ_codes_array);    %Avg. minimum no. of bits per symbol
bits_occ_codes = bits_occ_codes_persymbol*length(occ_codes_array);    %Total no. of bits for all symbols
disp(['Total entropy bits for occupancy codes: ' num2str(bits_occ_codes) ' (' num2str(bits_occ_codes_persymbol) ' bits per symbol)']);
disp(['Occupancy codes bpv (bits per input voxel): ' num2str(bits_occ_codes/size(ptcloud, 1))]);
disp(' ');

%Do the below for comparison only, of the number of bits required for
%unpruned vs pruned occupancy codes
if prune_flag == 1
    %Get the unpruned occupancy codes that would be transmitted to the 
    %decoder
    occupancy_codes_wout_pruning = cell(b, 1);
    for i = 1:(max_lvl - 1)
        occupancy_codes_wout_pruning{i} = myOT.OccupancyCode{i};
    end
    %Concatenate all of the occupancy codes (decimal values) at all octree 
    %levels from the root to (max_lvl - 1), into one long array
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
    disp(['Total entropy bits for unpruned occupancy codes: ' num2str(bits_unpruned_occ_codes) ' (' num2str(bits_unpruned_occ_codes_persymbol) ' bits per symbol)']);
    disp(['Unpruned occupancy codes bpv (bits per input voxel): ' num2str(bits_unpruned_occ_codes/size(ptcloud, 1))]);
    disp(' ');
    disp('*************************************');
    disp(' ');
end

%---- Post-Pruning Array (only if pruning has been done) ----

if prune_flag == 1    
    %Get the (parts of the) post_pruning_array that will be transmitted to the 
    %decoder
    post_pruning_array_forDec = cell(b, 1);
    for i = 1:(max_lvl - 1)
        post_pruning_array_forDec{i} = post_pruning_array{i};
    end
    %Concatenate all of the bits from post_pruning_array_forDec, at all octree 
    %levels from the root to (max_lvl - 1), into one long array
    pp_cntr = 1;
    pp_array = [];
    for l = 1:(max_lvl - 1)
        pp_array(pp_cntr:(pp_cntr + numel(post_pruning_array_forDec{l}) - 1)) = post_pruning_array_forDec{l};
        pp_cntr = pp_cntr + numel(post_pruning_array_forDec{l});
    end
    %Plot a histogram of the bits inside pp_array
    figure;
    histogram(pp_array);
    title({'Histogram of Bits Indicating Leaf/Non-Leaf Octree Cell', ['from Level 1-' num2str(max_lvl - 1)]});
    %Save the above histogram as a MATLAB figure and as a PDF image in our
    %network directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_pp_bits_lvl1-' num2str(max_lvl - 1)]);
    print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_pp_bits_lvl1-' num2str(max_lvl - 1)], '-dpdf');
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

%Get the quantized reconstructed control points (at start_lvl) that will be 
%transmitted to the decoder
rec_ctrlpts_forDec = reconstructed_control_points{start_lvl};
%Plot a histogram of the quantized control points inside rec_ctrlpts_forDec
figure;
histogram(rec_ctrlpts_forDec);
title(['Histogram of Quantized Control Points at Level ' num2str(start_lvl)]);
%Save the above histogram as a MATLAB figure and as a PDF image in our
%network directory (NB: The '-bestfit' option maximizes the size of the 
%figure to fill the page, but preserves the aspect ratio of the figure. 
%The figure might not fill the entire page. This option leaves a 
%minimum page margin of .25 inches).
savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_ctrlpts_lvl' num2str(start_lvl)]);
print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_ctrlpts_lvl' num2str(start_lvl)], '-dpdf');
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
%Plot a histogram of the quantized wavelet coefficients inside 
%wavelet_cfs_array
figure;
histogram(wavelet_cfs_array);
title(['Histogram of Quantized Wavelet Coefficients from Level ' num2str((start_lvl + 1)) '-' num2str(max_lvl)]);
%Save the above histogram as a MATLAB figure and as a PDF image in our
%network directory (NB: The '-bestfit' option maximizes the size of the 
%figure to fill the page, but preserves the aspect ratio of the figure. 
%The figure might not fill the entire page. This option leaves a 
%minimum page margin of .25 inches).
savefig(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_wavcfs_lvl' num2str(start_lvl + 1) '-' num2str(max_lvl)]);
print('-bestfit', ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\transmitted_histogram_quant_wavcfs_lvl' num2str(start_lvl + 1) '-' num2str(max_lvl)], '-dpdf');
%Compute the number of bits required for these wavelet coefficients
bits_wavelet_cfs_persymbol = entropy_calc(wavelet_cfs_array);  %Avg. minimum no. of bits per symbol
bits_wavelet_cfs = bits_wavelet_cfs_persymbol*length(wavelet_cfs_array);    %Total no. of bits for all symbols
disp(['Total entropy bits for wavelet coefficients: ' num2str(bits_wavelet_cfs) ' (' num2str(bits_wavelet_cfs_persymbol) ' bits per symbol)']);
disp(['Wavelet coefficients bpv (bits per input voxel): ' num2str(bits_wavelet_cfs/size(ptcloud, 1))]);
disp(' ');

%Do the below for comparison only, of the number of bits required for
%unpruned vs pruned wavelet coefficients
if prune_flag == 1
    %Get the unpruned quantized wavelet coefficients that would be 
    %transmitted to the decoder
    wavelet_coeffs_wout_pruning = cell(b, 1);
    for j = (start_lvl + 1):max_lvl
        wavelet_coeffs_wout_pruning{j} = wavelet_coeffs{j};
    end
    %Concatenate all of the wavelet coefficients at all octree levels from 
    %(start_lvl + 1) to max_lvl, into one long array
    wcf_cntr_wout_pruning = 1;
    wcf_array_wout_pruning = [];
    for l = (start_lvl + 1):max_lvl
        wcf_array_wout_pruning(wcf_cntr_wout_pruning:(wcf_cntr_wout_pruning + numel(wavelet_coeffs_wout_pruning{l}) - 1)) = wavelet_coeffs_wout_pruning{l};
        wcf_cntr_wout_pruning = wcf_cntr_wout_pruning + numel(wavelet_coeffs_wout_pruning{l});
    end
    %For comparison only, compute the number of bits that would be required 
    %for the unpruned wavelet coefficients
    disp('******** FOR COMPARISON ONLY ********');
    disp(' ');
    bits_unpruned_wcf_persymbol = entropy_calc(wcf_array_wout_pruning);    %Avg. minimum no. of bits per symbol
    bits_unpruned_wcf = bits_unpruned_wcf_persymbol*length(wcf_array_wout_pruning);    %Total no. of bits for all symbols
    disp(['Total entropy bits for unpruned wavelet coefficients: ' num2str(bits_unpruned_wcf) ' (' num2str(bits_unpruned_wcf_persymbol) ' bits per symbol)']);
    disp(['Unpruned wavelet coefficients bpv (bits per input voxel): ' num2str(bits_unpruned_wcf/size(ptcloud, 1))]);
    disp(' ');
    disp('*************************************');
    disp(' ');
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
disp(' ');
disp('------------------- ENCODER FINISHED -----------------------');
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










