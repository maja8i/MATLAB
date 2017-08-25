function [wavelet_coeffs, reconstructed_control_points] = wavelet_analysis(myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize)

%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) across all 
%octree blocks and levels, starting from start_lvl and going up to one
%level before the leaves
wavelet_coeffs = cell((b + 1), 1);  %These will be quantized coefficients
%Initialize a cell array to store the reconstructed signal at each vertex
%of each octree cell at every level from start_lvl to one level before the
%leaves
reconstructed_control_points = cell(size(control_points, 1), 1);

%Quantize all the control points at octree level start_lvl
for cpt = 1:size(control_points{start_lvl}, 1)
    control_points{start_lvl}(cpt, :) = quantize_uniform_scalar(control_points{start_lvl}(cpt, :), q_stepsize);
end
%Add these control points to reconstructed_control_points
reconstructed_control_points{start_lvl} = control_points{start_lvl};

%For each octree level, starting from start_lvl ...
for lvl = start_lvl:(max_lvl - 1)
    disp(['Computing wavelet coefficients between octree levels ' num2str(lvl) ' and ' num2str(lvl + 1) ' ...']);
    disp('------------------------------------------------------------');
    tic;
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Extract the current cell's 8 corner coordinates. This cell will 
        %represent our parent cell at the current level, and we will 
        %compute wavelet coefficients for its child cells.
        parent_corner_coords = corner_coords{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);  
        %Get the control points on all 8 corners of the current parent cell
        parent_ctrlpts_inds = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8));
        parent_control_points = control_points{lvl}(parent_ctrlpts_inds);
        
        %Find the midpoint coordinates of all the 12 edges of the current
        %parent block (we know in advance how the vertices are connected)
        parent_edge_midpoints = [(parent_corner_coords(1, :) + parent_corner_coords(2, :))/2;
            (parent_corner_coords(2, :) + parent_corner_coords(3, :))/2;
            (parent_corner_coords(3, :) + parent_corner_coords(4, :))/2;
            (parent_corner_coords(4, :) + parent_corner_coords(1, :))/2;
            (parent_corner_coords(1, :) + parent_corner_coords(5, :))/2;
            (parent_corner_coords(2, :) + parent_corner_coords(6, :))/2;
            (parent_corner_coords(3, :) + parent_corner_coords(7, :))/2;
            (parent_corner_coords(4, :) + parent_corner_coords(8, :))/2;
            (parent_corner_coords(5, :) + parent_corner_coords(6, :))/2;
            (parent_corner_coords(6, :) + parent_corner_coords(7, :))/2;
            (parent_corner_coords(7, :) + parent_corner_coords(8, :))/2;
            (parent_corner_coords(8, :) + parent_corner_coords(5, :))/2];
        
        %Find the average parent control points on the edge midpoints
        avg_edge_ctrlpts = [(parent_control_points(1) + parent_control_points(2))/2;
            (parent_control_points(2) + parent_control_points(3))/2;
            (parent_control_points(3) + parent_control_points(4))/2;
            (parent_control_points(4) + parent_control_points(1))/2;
            (parent_control_points(1) + parent_control_points(5))/2;
            (parent_control_points(2) + parent_control_points(6))/2;
            (parent_control_points(3) + parent_control_points(7))/2;
            (parent_control_points(4) + parent_control_points(8))/2;
            (parent_control_points(5) + parent_control_points(6))/2;
            (parent_control_points(6) + parent_control_points(7))/2;
            (parent_control_points(7) + parent_control_points(8))/2;
            (parent_control_points(8) + parent_control_points(5))/2];

        %Find the midpoint coordinates of all the 6 faces of the current 
        %parent block
        parent_face_midpoints = [mean([parent_corner_coords(1, :); parent_corner_coords(2, :); parent_corner_coords(5, :); parent_corner_coords(6, :)]);
            mean([parent_corner_coords(2, :); parent_corner_coords(3, :); parent_corner_coords(6, :); parent_corner_coords(7, :)]);
            mean([parent_corner_coords(3, :); parent_corner_coords(4, :); parent_corner_coords(7, :); parent_corner_coords(8, :)]);
            mean([parent_corner_coords(1, :); parent_corner_coords(4, :); parent_corner_coords(5, :); parent_corner_coords(8, :)]);
            mean([parent_corner_coords(1, :); parent_corner_coords(2, :); parent_corner_coords(3, :); parent_corner_coords(4, :)]);
            mean([parent_corner_coords(5, :); parent_corner_coords(6, :); parent_corner_coords(7, :); parent_corner_coords(8, :)])];
        
        %Find the average parent control points on the face midpoints
        avg_face_ctrlpts = [mean([parent_control_points(1); parent_control_points(2); parent_control_points(5); parent_control_points(6)]);
            mean([parent_control_points(2); parent_control_points(3); parent_control_points(6); parent_control_points(7)]);
            mean([parent_control_points(3); parent_control_points(4); parent_control_points(7); parent_control_points(8)]);
            mean([parent_control_points(1); parent_control_points(4); parent_control_points(5); parent_control_points(8)]);
            mean([parent_control_points(1); parent_control_points(2); parent_control_points(3); parent_control_points(4)]);
            mean([parent_control_points(5); parent_control_points(6); parent_control_points(7); parent_control_points(8)])];

        %Get the pointers to all the child cells of the current parent cell
        %(occ_cell)
        child_ptrs = myOT.FirstChildPtr{lvl}(occ_cell):(myOT.FirstChildPtr{lvl}(occ_cell) + uint32(myOT.ChildCount{lvl}(occ_cell)) - 1);
        %Get the corner coordinates of all the children of the current
        %parent cell
        cnr_coords_inds = ((8*child_ptrs(1) - 7):(8*child_ptrs(end)))'; %Indices for corner_coords{lvl + 1} for the current children
        child_corner_coords = corner_coords{lvl + 1}(cnr_coords_inds, :);
        %Get the control points of all the child corners
        child_ctrlpts_inds = ctrl_pts_pointers{lvl + 1}(cnr_coords_inds);
        child_control_points = control_points{lvl + 1}(child_ctrlpts_inds);
        
        %Initialize an array that will store, for each child corner, the 
        %average of the control points on its corresponding parent corners
        averages = zeros(size(child_corner_coords, 1), 1);
        
        %For each corner of the current children, check where the corner is
        %relative to the parent block:   
        
        %On a parent corner (i.e., same as a parent corner)
        on_parent_child_inds = find(ismember(child_corner_coords, parent_corner_coords, 'rows') > 0);
        %In this case, the signal on each of these corners is a low-pass 
        %coefficient (not a wavelet coefficient) and has already been 
        %reconstructed as it is equal to its corresponding parent control 
        %point. So, do nothing but place the value of the parent control 
        %point in the corresponding location in the "averages" array - this 
        %will ensure that the wavelet coefficient for the corresponding 
        %child will be equal to 0 and the reconstructed control point for 
        %the child will be the same as the parent control point.
        if ~isempty(on_parent_child_inds)
            p_inds = zeros(length(on_parent_child_inds), 1);
            for i = 1:length(on_parent_child_inds) 
                %Find which parent corner corresponds to which child corner
                %represented in on_parent_child_inds
                p_inds(i) = find(ismember(parent_corner_coords, child_corner_coords(on_parent_child_inds(i), :), 'rows') > 0);
            end
            averages(on_parent_child_inds) = parent_control_points(p_inds);
        end
   
        %On a parent edge
        on_edge_child_inds = find(ismember(child_corner_coords, parent_edge_midpoints, 'rows') > 0);
        %In this case, compute the average of the Bezier control points 
        %found on the 2 corners of the corresponding parent edge
        if ~isempty(on_edge_child_inds)
            p_inds = zeros(length(on_edge_child_inds), 1);
            for i = 1:length(on_edge_child_inds) 
                %Find which parent edge midpoint corresponds to which child 
                %corner represented in on_edge_child_inds
                p_inds(i) = find(ismember(parent_edge_midpoints, child_corner_coords(on_edge_child_inds(i), :), 'rows') > 0);
            end
            averages(on_edge_child_inds) = avg_edge_ctrlpts(p_inds);
        end
        
        %On a parent face
        on_face_child_inds = find(ismember(child_corner_coords, parent_face_midpoints, 'rows') > 0);
        %In this case, compute the average of the Bezier control points 
        %found on the 4 corners of the corresponding parent face
        if ~isempty(on_face_child_inds)
            p_inds = zeros(length(on_face_child_inds), 1);
            for i = 1:length(on_face_child_inds) 
                %Find which parent face midpoint corresponds to which child 
                %corner represented in on_face_child_inds
                p_inds(i) = find(ismember(parent_face_midpoints, child_corner_coords(on_face_child_inds(i), :), 'rows') > 0);
            end
            averages(on_face_child_inds) = avg_face_ctrlpts(p_inds);
        end
        
        %In the centre of the parent block
        in_centre_child_inds = find(ismember(child_corner_coords, mean(parent_corner_coords), 'rows'));
        %In this case, compute the average of the Bezier control points 
        %found on all 8 corners of the current parent block
        if ~isempty(in_centre_child_inds)
            averages(in_centre_child_inds) = mean(parent_control_points);
        end
        
        %For all the child corners of the current parent block, subtract 
        %their corresponding average parent signal (control point) computed
        %above, from the child control point. The result will be the high-
        %pass transform (wavelet) coefficient of the child corner.
        wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = child_control_points - averages;
        %Quantize the wavelet coefficients computed above
        wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1), q_stepsize);
        %Add the quantized wavelet coefficients to the corresponding values
        %in "averages", to obtain the reconstructed signal (control point) 
        %at each corresponding child corner
        reconstructed_control_points{lvl + 1}(cnr_coords_inds, 1) = averages + wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1);
        
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
        
        %disp(['Finished occ_cell ' num2str(occ_cell) '/' num2str(myOT.NodeCount(lvl))]);
    end %End occ_cell
    
    %Keep only the wavelet coefficients and reconstructed control points
    %for the UNIQUE child corners at level lvl + 1, and discard the rest
    wavelet_coeffs{lvl + 1} = wavelet_coeffs{lvl + 1}(unique(ctrl_pts_pointers{lvl + 1}, 'stable'));
    reconstructed_control_points{lvl + 1} = reconstructed_control_points{lvl + 1}(unique(ctrl_pts_pointers{lvl + 1}, 'stable'));
    
    wcfs_time = toc;
    disp(['Time taken to compute wavelet coefficients and reconstructed control points for level ' num2str(lvl + 1) ': ' num2str(wcfs_time) ' seconds']);
    disp('------------------------------------------------------------');
end %End lvl