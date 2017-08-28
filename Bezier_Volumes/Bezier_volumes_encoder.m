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
[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted] = construct_octree(ptcloud_file, b);
%[myOT, mortonCodes_sorted, xyz_sorted, normals_sorted] = construct_octree(ptcloud_file, b);

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
[~, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_normal_averages, occupied_voxel_centroids, occupied_voxel_centroid_averages] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, normals_sorted, centroids_sorted);
%[~, occupied_voxel_coords, occupied_voxel_normals] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, normals_sorted);

disp(' ');
disp('------------ Computing Bezier Control Points ---------------');
disp(' ');

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

%Set a threshold, above which control points will be computed as average
%dot products for each unique corner.  The difference vectors in this case
%will be measured as distances between corner coordinates and either voxel 
%centroids (at the voxel level) or the average centroid of all the occupied 
%voxels associated with a given octree block (at levels other than the 
%voxel level). To obtain the dot product for each corner, these difference 
%vectors will be projected onto the unit-norm normal vector of the 
%corresponding voxel (at the voxel level) or the unit-norm average normal 
%of all the occupied voxels associated with a given octree block (at levels
%other than the voxel level). The control point for each UNIQUE corner will
%then be computed as the average of all the dot products computed for the
%same corner (for different voxels or blocks that share that corner). For
%octree levels below this threshold, the difference vector for each unique 
%corner will be computed as the Euclidean distance between this corner and
%the nearest voxel's midpoint, and the absolute value of this minimum 
%distance will be the absolute value of the control point for the 
%corresponding corner. The sign of each control point in this case will be
%obtained as the sign of the dot product between the corresponding 
%difference vector and the normal vector of the corresponding nearest voxel
%found.
thresh = 9; %Root is level 1 in MATLAB, so thresh corresponds to level 8 when root is 0
    
%tic;
%for lvl = 1:(b + 1)
for lvl = 1:max_lvl
    disp(['Computing Bezier control points for octree level ' num2str(lvl) ':']); 
    tic;
    if lvl < thresh
        %For each unique corner coordinate, find the nearest occupied voxel 
        %in any of the occupied octree cells at the current level that 
        %share this corner
        for c = 1:size(unique_coords{lvl}, 1)
            %Find out which occupied octree cells share this corner
            shared_cells{lvl, c} = ceil(find((c - ctrl_pts_pointers{lvl}) == 0)./8);    %Faster than using "ismember"                
            %Get the (x, y, z) coordinates of all the occupied voxels in the 
            %cells found above, and the normals for each of these voxels
            current_occupied_voxel_coords = occupied_voxel_coords{lvl, shared_cells{lvl, c}(1)};    %Coords. are voxel midpoints 
            current_occupied_voxel_normals = occupied_voxel_normals{lvl, shared_cells{lvl, c}(1)};
            if length(shared_cells{lvl, c}) > 1
                %Concatenate the position and normal coordinates from the cells 
                %that share the current corner (concatenate along rows, so that
                %we maintain the 3-column structure)
                current_occupied_voxel_coords = cat(1, current_occupied_voxel_coords, occupied_voxel_coords{lvl, shared_cells{lvl, c}(2:length(shared_cells{lvl, c}))}); %Coords. are voxel midpoints
                current_occupied_voxel_normals = cat(1, current_occupied_voxel_normals, occupied_voxel_normals{lvl, shared_cells{lvl, c}(2:length(shared_cells{lvl, c}))});
            end
            %Find the nearest neighbour to the current unique corner coordinate 
            diff = unique_coords{lvl}(c, :) - current_occupied_voxel_coords; 
            %euclid_dists = abs(arrayfun(@(idx) norm(diff(idx, :)), 1:size(diff, 1)));
            %Don't need "abs" below, because the diff values are squared
            euclid_dists = sqrt(sum(diff.^2, 2));   %This operation is faster than "arrayfun", above
            [min_euclid_dist{lvl}(c), nearest_voxel_ind] = min(euclid_dists);
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
            if (lvl > 1)&&(lvl < b + 1)
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
            %disp(['Finished corner ' num2str(c) '/' num2str(size(unique_coords{lvl}, 1))]);
        end %End c (current unique coordinate)      
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
    elseif lvl >= thresh   
        %If we are at the voxel level ...
        if lvl == b + 1
            %Turn the occupied_voxel_centroids cell array at this level, 
            %into a matrix: we only have one centroid per voxel anyway, so
            %collect all of the voxels' centroids into this one matrix
            occ_vox_centroids = cat(1, occupied_voxel_centroids{lvl, 1}, occupied_voxel_centroids{lvl, 2:end});
            %Compute the difference between each corner at this level and
            %the centroid vector of the corresponding voxel (some corners 
            %will be used by more than one voxel, so the calculation will 
            %be repeated for each shared voxel)
            diff = corner_coords{lvl} - occ_vox_centroids(repmat(1:size(occ_vox_centroids, 1), 8, 1), :);   %Repeat each voxel centroid 8 times
            %Turn the occupied_voxel_normals cell array at this level, into 
            %a matrix: we only have one normal vector per voxel anyway, so
            %collect all of the voxels' normals into this one matrix
            occ_vox_normals = cat(1, occupied_voxel_normals{lvl, 1}, occupied_voxel_normals{lvl, 2:end});
            %Compute the dot product between each diff vector and the
            %corresponding voxel's unit-norm normal vector (this is the
            %scalar projection of each difference vector onto the 
            %corresponding voxel's normal vector)
            unit_norm_normals = (occ_vox_normals(repmat(1:size(occ_vox_normals, 1), 8, 1), :)./(sqrt(sum((occ_vox_normals(repmat(1:size(occ_vox_normals, 1), 8, 1), :)).^2, 2))));  %Repeat each voxel normal 8 times
            dot_prods = dot(diff, unit_norm_normals, 2);
        %If we are not at the voxel level ...
        else
            %Turn the occupied_voxel_centroid_averages cell array at this 
            %level, into a matrix: we have one centroid average per
            %occupied block, so collect all of these averages into this one
            %matrix
            occ_vox_centroid_avgs = cat(1, occupied_voxel_centroid_averages{lvl, 1}, occupied_voxel_centroid_averages{lvl, 2:myOT.NodeCount(lvl)});
            %Compute the difference between each corner at this level and
            %the average centroid vector of the corresponding octree block
            %(some corners will be used by more than one block, so the 
            %calculation will be repeated for each shared block)
            diff = corner_coords{lvl} - occ_vox_centroid_avgs(repmat(1:size(occ_vox_centroid_avgs, 1), 8, 1), :);   %Repeat each average centroid 8 times
            %Turn the occupied_voxel_normal_averages cell array at this 
            %level, into a matrix: we have one normal vector average per
            %occupied block, so collect all of these averages into this one
            %matrix
            occ_vox_normal_avgs = cat(1, occupied_voxel_normal_averages{lvl, 1}, occupied_voxel_normal_averages{lvl, 2:myOT.NodeCount(lvl)});
            %Compute the dot product between each diff vector and the
            %corresponding octree block's unit-norm average normal vector 
            %(this is the scalar projection of each difference vector onto 
            %the corresponding block's average normal vector)
            unit_norm_normals = (occ_vox_normal_avgs(repmat(1:size(occ_vox_normal_avgs, 1), 8, 1), :)./(sqrt(sum((occ_vox_normal_avgs(repmat(1:size(occ_vox_normal_avgs, 1), 8, 1), :)).^2, 2))));  %Repeat each average normal 8 times
            dot_prods = dot(diff, unit_norm_normals, 2);
        end
        %Create an M x N matrix A, where M = size(unique_coords{lvl}, 1) 
        %and each row represents one unique coordinate, and
        %N = size(corner_coords{lvl}, 1) so that each column represents one
        %corner coordinate out of all the corner coordinates at the current
        %octree level (not just the unique coordinates). A(n, m) = 1 if
        %corner_coords{lvl}(n, :) == unique_coords{lvl}(m, :), and = 0 
        %otherwise.
        M = size(unique_coords{lvl}, 1); %Number of unique corners at the current octree level
        N = size(corner_coords{lvl}, 1); %Total number of corners at the current octree level, including repetitions from all blocks
        A = sparse(ctrl_pts_pointers{lvl}, [1:N]', ones(N,1), M, N);
        %Compute the control point for each UNIQUE corner at the current
        %level as the average of all the dot products computed for the
        %same corner (i.e., the average of all the dot_prods associated
        %with the same unique corner)
        uniqueCount = full(sum(A, 2));
        uniqueTotal = A*dot_prods;
        control_points{lvl} = uniqueTotal./(repmat(uniqueCount, 1, size(uniqueTotal, 2))); %Valid only when uniqueCount > 0, else NaN        
    end %End check if lvl < thresh
    cp_time = toc;
    disp(' ');
    disp(['Time taken to compute control points at level ' num2str(lvl) ': ' num2str(cp_time) ' seconds']);
    disp('------------------------------------------------------------');
end %End lvl
% ctrlpts_time = toc;
% disp(' ');
% disp('************************************************************');
% disp(['Time taken to compute all control points: ' num2str(ctrlpts_time) ' seconds']);
% disp('************************************************************');

% %For debugging purposes, print out how many voxels, if any, end up having 
% %all 8 of their control points with the same sign (+/-) ...
% [~, ptcloud, ~] = plyRead(ptcloud_file);
% %Accumulate all the control points (not just the unique ones) into one cell
% %array
% all_ctrlpts = get_all_ctrlpts(control_points, ctrl_pts_pointers, start_lvl, max_lvl); 
% same_sign_cntr = 0;
% vox = 0;
% same_sign_voxels = [];
% disp(' ');
% for i = 1:8:(length(all_ctrlpts{max_lvl}) - (max_lvl - 1)) 
%     vox = vox + 1;
%     %Get all 8 control points for the corners of the current voxel
%     current_ctrlpts = all_ctrlpts{max_lvl}(i:(i + max_lvl - 1));
%     %Check the signs of current_ctrlpts (assume here that none of the
%     %control points have a value of 0)
%     if (sum(sign(current_ctrlpts)) == length(current_ctrlpts))||(sum(sign(current_ctrlpts)) == -length(current_ctrlpts))
%         same_sign_cntr = same_sign_cntr + 1;
% %         disp(' ');
%         disp(['Voxel ' num2str(vox) ' (' num2str(ptcloud(vox, 1:3)) ') has all control points with the same sign: ']);
% %         disp(' ');
%         %Store this voxel, for debugging purposes
%         same_sign_voxels((size(same_sign_voxels, 1) + 1), 1:3) = ptcloud(vox, 1:3);
% %         %Get the corner coordinates for each corner of this voxel
% %         vox_corners = corner_coords{max_lvl}((i:(i + max_lvl - 1)), 1:3);
% %         disp('Voxel corner coordinates:');
% %         disp(num2str(vox_corners));
% %         disp(' ');
% %         %Get the coordinates (either midpoint or centroid, whichever one of
% %         %these happened to be used in the control point computation) of the
% %         %nearest voxel found for each corner of the current voxel
% %         curr_nearest_vox = nearest_voxels{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
% %         disp('Nearest voxel found for each corner of the current voxel: ');
% %         disp(num2str(curr_nearest_vox));
% %         disp(' ');
% %         %Get the normal vector for the current voxel, from the input
% %         %voxelized point cloud
% %         disp('Voxel normal:');
% %         disp(normals_sorted(vox, :));
% %         %Get the normal vector for each nearest voxel found above
% %         curr_normal_vec = normal_nearest_vox{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
% %         disp('Nearest voxel normal for each corner: ');
% %         disp(num2str(curr_normal_vec));
% %         disp(' ');
% %         %Get the difference vector for each of the 8 corners of this voxel
% %         %(difference vector between each corner and the chosen voxel)
% %         vox_diff_vecs = difference_vectors{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)), 1:3);
% %         disp('Voxel corner difference vectors:');
% %         disp(num2str(vox_diff_vecs));
% %         disp(' ');
% %         %Get the dot product between the difference vector of each corner
% %         %of this voxel and the normal vector of the chosen nearest voxel
% %         %(the chosen voxel may not be the current voxel itself, as
% %         %neighbouring voxels are the same distance away)
% %         vox_dp = dot_products{max_lvl}(ctrl_pts_pointers{max_lvl}(i:(i + max_lvl - 1)));
% %         disp('Voxel corner dot products:');
% %         disp(num2str(vox_dp));
% %         disp(' ');
%     end
% end
% disp(['TOTAL number of voxels with all control points having the same sign: ' num2str(same_sign_cntr) '/' num2str(length(all_ctrlpts{max_lvl})/8) ' (' num2str((same_sign_cntr/(length(all_ctrlpts{max_lvl})/8))*100) '%)']);
% %Plot voxels that have the same control point signs
% if same_sign_cntr > 0
%     figure;
%     scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, 'filled', 'MarkerFaceColor', 'b');
%     hold on;
%     scatter3(same_sign_voxels(:, 1), same_sign_voxels(:, 2), same_sign_voxels(:, 3), 5, 'filled', 'MarkerFaceColor', 'm');
%     axis equal; axis off;
%     title('Voxels with Same-Sign Control Points at Encoder');
%     legend('Voxels with Different-Sign Control Points', 'Voxels with Same-Sign Control Points', 'Location', 'best');
% end

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

[wavelet_coeffs, reconstructed_control_points] = wavelet_analysis(myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize);

%Plot the distribution of reconstructed control points (low-pass
%coefficients) at different octree levels
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

%------------------------------- Pruning  --------------------------------%

%Prune octree and wavelet coefficient tree, if pruning is required
if prune_flag == 1
    disp(' ');
    disp('-------------------- Pruning Octree ------------------------');
    disp(' ');

    [pruned_occupancy_codes, post_pruning_array, toprune, toprune2] = prune_octree(myOT, all_zero_wav_cfs, start_lvl, max_lvl, b);

    disp(' ');
    disp('------------ Pruning Wavelet Coefficient Tree --------------');
    disp(' ');

    pruned_wavelet_coeffs = prune_wavelet_coeff_tree(wavelet_coeffs, toprune, toprune2, shared_cells, ctrl_pts_pointers, myOT, start_lvl, max_lvl, b);
end
 
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
    if max_lvl == b + 1
        for i = 1:(max_lvl - 1) %Assumes max_lvl = b + 1
            post_pruning_array_forDec{i} = post_pruning_array{i};
        end
    elseif max_lvl < b + 1
        for i = 1:max_lvl
            post_pruning_array_forDec{i} = post_pruning_array{i};
        end
    end
    %Concatenate all of the bits from post_pruning_array_forDec, at all octree 
    %levels from the root to (max_lvl - 1), into one long array
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










