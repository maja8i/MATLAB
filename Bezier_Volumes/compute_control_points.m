function [control_points, nearest_voxels, normal_nearest_vox, min_euclid_dist, difference_vectors, dot_products, thresh] = compute_control_points(debug_flag, myOT, corner_coords, unique_coords, ctrl_pts_pointers, occupied_voxel_coords, occupied_voxel_normals, occupied_voxel_centroids, occupied_voxel_normal_averages, occupied_voxel_centroid_averages, b, max_lvl)

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
%Initialize a cell array to store the Bezier control points (signed
%distances) associated with the unique corner coordinates at each octree 
%level
control_points = cell((b + 1), 1);

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
thresh = 9; %Root is level 1 in MATLAB
    
start_ctrlpts_time = tic;
%for lvl = 1:(b + 1)
for lvl = 1:max_lvl
    if debug_flag == 1
        disp(['Computing Bezier control points for octree level ' num2str(lvl) ':']); 
        start_cp_time_lvl = tic;
    end
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
        control_points{lvl} = uniqueTotal./uniqueCount; %Valid only when uniqueCount > 0, else NaN        
    end %End check if lvl < thresh
    if debug_flag == 1
        cp_time_lvl = toc(start_cp_time_lvl);
        disp(' ');
        disp(['Time taken to compute control points at level ' num2str(lvl) ': ' num2str(cp_time_lvl) ' seconds']);
        disp('------------------------------------------------------------');
    end
end %End lvl
ctrlpts_time = toc(start_ctrlpts_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all control points: ' num2str(ctrlpts_time) ' seconds']);
disp('************************************************************');