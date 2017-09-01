function reconstruction_decoder = reconstruct_control_points_decoder(rec_ctrlpts_forDec, wavelet_coeffs_forDec, OccupancyCode, FirstChildPtr, ChildCount, corner_coords_decoder, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, prune_flag, reconstructed_control_points)

%Assign inputs to new variables
rec_ctrlpts_start_lvl = rec_ctrlpts_forDec;
wavelet_coeffs = wavelet_coeffs_forDec;

%Initialize a cell array to store the reconstructed signal (control point)
%at each vertex of each octree cell at every level from start_lvl to the 
%leaves
reconstruction_decoder = cell((b + 1), 1);

%Dequantize the reconstructed control points for start_lvl, which were 
%received from the encoder, then place them in the corresponding cell of 
%reconstruction_decoder. Note that when q_stepsize = 1, the dequantized
%values will just be equal to the quantized values.
reconstruction_decoder{start_lvl} = dequantize_uniform_scalar(rec_ctrlpts_start_lvl, q_stepsize);

%Start from the start_lvl and work up to the leaves ...
%tic;
%for lvl = start_lvl:b
for lvl = start_lvl:(max_lvl - 1)
    disp(['Reconstructing control points at octree level ' num2str(lvl + 1) ' ...']);
    tic;
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %Dequantize the wavelet coefficients for all the corners (not just the
    %unique ones) of all the occupied children (at level lvl + 1) of all 
    %the occupied octree cells at the current level (lvl) - that is, get
    %the dequantized wavelet coefficients for all the occupied octree cells
    %at level lvl + 1
    all_wavelet_coeffs = dequantize_uniform_scalar(wavelet_coeffs{lvl + 1}(ctrl_pts_pointers{lvl + 1}), q_stepsize);

    %For each occupied octree cell at the current level ...
    for occ_cell = 1:numel(OccupancyCode{lvl})
        %Extract the cell's 8 corner coordinates. This cell will represent
        %our parent cell at the current level, since we will consider its
        %children for the control point (signal) reconstruction.
        parent_corner_coords = corner_coords_decoder{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);    
        %Find the average parent x, y, z coordinates
        parent_avg_coords = sum(parent_corner_coords)./8;
        %Get the reconstructed control points on all 8 corners of the 
        %current parent cell
        parent_ctrlpts_inds = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8));
        parent_control_points = reconstruction_decoder{lvl}(parent_ctrlpts_inds);        
        %Find the midpoint coordinates of all the 12 edges of the current
        %parent block (we know in advance how the vertices are connected)
        parent_edge_midpoints = [(parent_corner_coords(1, :) + parent_corner_coords(2, :));
            (parent_corner_coords(2, :) + parent_corner_coords(3, :));
            (parent_corner_coords(3, :) + parent_corner_coords(4, :));
            (parent_corner_coords(4, :) + parent_corner_coords(1, :));
            (parent_corner_coords(1, :) + parent_corner_coords(5, :));
            (parent_corner_coords(2, :) + parent_corner_coords(6, :));
            (parent_corner_coords(3, :) + parent_corner_coords(7, :));
            (parent_corner_coords(4, :) + parent_corner_coords(8, :));
            (parent_corner_coords(5, :) + parent_corner_coords(6, :));
            (parent_corner_coords(6, :) + parent_corner_coords(7, :));
            (parent_corner_coords(7, :) + parent_corner_coords(8, :));
            (parent_corner_coords(8, :) + parent_corner_coords(5, :))]./2;
        %Find the average parent control points on the edge midpoints
        avg_edge_ctrlpts = [(parent_control_points(1) + parent_control_points(2));
            (parent_control_points(2) + parent_control_points(3));
            (parent_control_points(3) + parent_control_points(4));
            (parent_control_points(4) + parent_control_points(1));
            (parent_control_points(1) + parent_control_points(5));
            (parent_control_points(2) + parent_control_points(6));
            (parent_control_points(3) + parent_control_points(7));
            (parent_control_points(4) + parent_control_points(8));
            (parent_control_points(5) + parent_control_points(6));
            (parent_control_points(6) + parent_control_points(7));
            (parent_control_points(7) + parent_control_points(8));
            (parent_control_points(8) + parent_control_points(5))]./2;        
        %Find the midpoint coordinates of all the 6 faces of the current 
        %parent block
        parent_face_midpoints = [(parent_corner_coords(1, :) + parent_corner_coords(2, :) + parent_corner_coords(5, :) + parent_corner_coords(6, :));
            (parent_corner_coords(2, :) + parent_corner_coords(3, :) + parent_corner_coords(6, :) + parent_corner_coords(7, :));
            (parent_corner_coords(3, :) + parent_corner_coords(4, :) + parent_corner_coords(7, :) + parent_corner_coords(8, :));
            (parent_corner_coords(1, :) + parent_corner_coords(4, :) + parent_corner_coords(5, :) + parent_corner_coords(8, :));
            (parent_corner_coords(1, :) + parent_corner_coords(2, :) + parent_corner_coords(3, :) + parent_corner_coords(4, :));
            (parent_corner_coords(5, :) + parent_corner_coords(6, :) + parent_corner_coords(7, :) + parent_corner_coords(8, :))]./4;
        %Find the average parent control points on the face midpoints
        avg_face_ctrlpts = [(parent_control_points(1) + parent_control_points(2) + parent_control_points(5) + parent_control_points(6));
            (parent_control_points(2) + parent_control_points(3) + parent_control_points(6) + parent_control_points(7));
            (parent_control_points(3) + parent_control_points(4) + parent_control_points(7) + parent_control_points(8));
            (parent_control_points(1) + parent_control_points(4) + parent_control_points(5) + parent_control_points(8));
            (parent_control_points(1) + parent_control_points(2) + parent_control_points(3) + parent_control_points(4));
            (parent_control_points(5) + parent_control_points(6) + parent_control_points(7) + parent_control_points(8))]./4;
        
        %Get the pointers to all the child cells of the current parent cell
        %(occ_cell)
        child_ptrs = FirstChildPtr{lvl}(occ_cell):(FirstChildPtr{lvl}(occ_cell) + uint32(ChildCount{lvl}(occ_cell)) - 1);
        %Get the corner coordinates of all the children of the current
        %parent cell
        cnr_coords_inds = ((8*child_ptrs(1) - 7):(8*child_ptrs(end)))'; %Indices for corner_coords{lvl + 1} for the current children
        child_corner_coords = corner_coords_decoder{lvl + 1}(cnr_coords_inds, :);
        
        %Initialize an array that will store, for each child corner, the 
        %average of the control points on its corresponding parent corners
        averages = zeros(size(child_corner_coords, 1), 1);
        
        %For each corner of the current children, check where the corner is
        %relative to the parent block:   
        
        %On a parent corner (i.e., same as a parent corner)
        sub = repmat(child_corner_coords, 8, 1) - parent_corner_coords(repmat(1:size(parent_corner_coords, 1), size(child_corner_coords, 1), 1), :);
        section_nbr = ceil((find(sum(abs(sub), 2) == 0)./size(child_corner_coords, 1)));
        on_parent_child_inds = find(sum(abs(sub), 2) == 0) - size(child_corner_coords, 1)*(section_nbr - 1);   
        %In this case, the signal on each of these corners is a low-pass 
        %coefficient (not a wavelet coefficient) and has already been 
        %reconstructed as it is equal to its corresponding parent control 
        %point. So, do nothing but place the value of the parent control 
        %point in the corresponding location in the "averages" array - this 
        %will ensure that the wavelet coefficient for the corresponding 
        %child will be equal to 0 and the reconstructed control point for 
        %the child will be the same as the parent control point.
        if ~isempty(on_parent_child_inds)
            averages(on_parent_child_inds) = parent_control_points(section_nbr);
        end
        
        %On a parent edge
        sub = repmat(child_corner_coords, size(parent_edge_midpoints, 1), 1) - parent_edge_midpoints(repmat(1:size(parent_edge_midpoints, 1), size(child_corner_coords, 1), 1), :);
        section_nbr = ceil((find(sum(abs(sub), 2) == 0)./size(child_corner_coords, 1)));
        on_edge_child_inds = find(sum(abs(sub), 2) == 0) - size(child_corner_coords, 1)*(section_nbr - 1);
        %In this case, compute the average of the Bezier control points 
        %found on the 2 corners of the corresponding parent edge
        if ~isempty(on_edge_child_inds)
            averages(on_edge_child_inds) = avg_edge_ctrlpts(section_nbr); 
        end
       
        %On a parent face
        sub = repmat(child_corner_coords, size(parent_face_midpoints, 1), 1) - parent_face_midpoints(repmat(1:size(parent_face_midpoints, 1), size(child_corner_coords, 1), 1), :);
        section_nbr = ceil((find(sum(abs(sub), 2) == 0)./size(child_corner_coords, 1)));
        on_face_child_inds = find(sum(abs(sub), 2) == 0) - size(child_corner_coords, 1)*(section_nbr - 1);
        %In this case, compute the average of the Bezier control points 
        %found on the 4 corners of the corresponding parent face
        if ~isempty(on_face_child_inds)
            averages(on_face_child_inds) = avg_face_ctrlpts(section_nbr); 
        end

        %In the centre of the parent block
        in_centre_child_inds = find(sum(abs(child_corner_coords - parent_avg_coords), 2) == 0);
        %In this case, compute the average of the Bezier control points 
        %found on all 8 corners of the current parent block
        if ~isempty(in_centre_child_inds)
            averages(in_centre_child_inds) = mean(parent_control_points);
        end
        
        %Get the dequantized wavelet coefficients for all the corners (not
        %just the unique ones) of all the children of the current parent 
        %cell
        child_wavelet_coeffs = all_wavelet_coeffs(cnr_coords_inds, 1);
        %Add the dequantized wavelet coefficients above, to the 
        %corresponding values in "averages", to obtain the reconstructed 
        %control point at each corresponding child corner
        reconstruction_decoder{lvl + 1}(cnr_coords_inds, 1) = child_wavelet_coeffs + averages;
        
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
    end %End occ_cell
    
    %Keep only the reconstructed control points for the UNIQUE child 
    %corners at level lvl + 1, and discard the rest
    [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers{lvl + 1}, 'stable');
    reconstruction_decoder{lvl + 1} = reconstruction_decoder{lvl + 1}(unique_ctrl_pts_pointers_inds, 1);

    cp_time = toc;
    disp(' ');
    disp(['Time taken to reconstruct control points at level ' num2str(lvl + 1) ': ' num2str(cp_time) ' seconds']);
    disp('------------------------------------------------------------');
    
    %BELOW IS FOR DEBUGGING PURPOSES ONLY:
    if prune_flag == 0
        %Check if the control points at this level have been correctly
        %reconstructed (i.e., if they are identical to the control points
        %at the encoder)
        test_ctrlpts = reconstructed_control_points{lvl + 1} - reconstruction_decoder{lvl + 1};
        test_ctrlpts_wrong = find(test_ctrlpts ~= 0);
        if isempty(test_ctrlpts_wrong)
            disp('All control points reconstructed correctly');
        else
            disp(['Number of incorrectly reconstructed control points: ' num2str(length(test_ctrlpts_wrong))]);
        end
    end
end %End lvl
% ctrlpt_recon_time = toc;
% disp(' ');
% disp('************************************************************');
% disp(['Time taken to reconstruct all control points: ' num2str(ctrlpt_recon_time) ' seconds']);
% disp('************************************************************');