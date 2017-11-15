function [wavelet_coeffs, reconstructed_control_points] = wavelet_analysis(debug_flag, myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, zero_threshold)

%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) across all 
%octree blocks and levels, starting from start_lvl and going up to one
%level before the leaves
wavelet_coeffs = cell((b + 1), 1);  %These will be quantized coefficients
%Initialize a cell array to store the reconstructed signal at each vertex
%of each octree cell at every level from start_lvl to one level before the
%leaves
reconstructed_control_points = cell(size(control_points, 1), 1);

start_wcfs_time = tic;

%Quantize and dequantize (reconstruct) all the control points at octree 
%level start_lvl. The dequantized values will be used for wavelet analysis
%and control point reconstruction, below, but the quantized values will be
%transmitted to the decoder.
quantized_cp = quantize_uniform_scalar(control_points{start_lvl}, q_stepsize);
reconstructed_control_points{start_lvl} = dequantize_uniform_scalar(quantized_cp, q_stepsize);
    
%For each octree level, starting from start_lvl ...
for lvl = start_lvl:(max_lvl - 1) 
    %profile on
    if debug_flag == 1
        disp(['Computing wavelet coefficients between octree levels ' num2str(lvl) ' and ' num2str(lvl + 1) ' ...']);
        disp('------------------------------------------------------------');
        start_wcfs_time_lvl = tic;
    end
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Extract the current cell's 8 corner coordinates. This cell will 
        %represent our parent cell at the current level, and we will 
        %compute wavelet coefficients for its child cells.
        parent_corner_coords = corner_coords{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);  
        %Find the average parent x, y, z coordinates
        parent_avg_coords = sum(parent_corner_coords)./8;
        %Get the (reconstructed) control points on all 8 corners of the current parent cell
        parent_ctrlpts_inds = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):(occ_cell*8));
        parent_control_points = reconstructed_control_points{lvl}(parent_ctrlpts_inds);
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
        
        %For all the child corners of the current parent block, subtract 
        %their corresponding average parent signal (control point) computed
        %above, from the child control point. The result will be the high-
        %pass transform (wavelet) coefficient of the child corner.
        wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = child_control_points - averages;
        %Quantize the wavelet coefficients computed above
        wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1), q_stepsize);
        %Check if any of the quantized wavelet coefficients are within +/-
        %zero_threshold of 0: if so, set the values of these quantized
        %coefficients to 0
        thresholded_wav_cfs = wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1);
        thresholded_wav_cfs(abs(thresholded_wav_cfs) <= zero_threshold) = 0;
        wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = thresholded_wav_cfs;
        %Dequantize the quantized wavelet coefficients and add them to the
        %corresponding values in "averages", to obtain the reconstructed 
        %signal (control point) at each corresponding child corner. Note
        %that the QUANTIZED versions of these wavelet coefficients will be
        %sent to the decoder.
        dequant_wavelet_coeffs = dequantize_uniform_scalar(wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1), q_stepsize);
        reconstructed_control_points{lvl + 1}(cnr_coords_inds, 1) = averages + dequant_wavelet_coeffs;
        
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
        
        %disp(['Finished occ_cell ' num2str(occ_cell) '/' num2str(myOT.NodeCount(lvl))]);
    end %End occ_cell
    
    %Keep only the wavelet coefficients and reconstructed control points
    %for the UNIQUE child corners at level lvl + 1, and discard the rest
    [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers{lvl + 1}, 'stable');
    wavelet_coeffs{lvl + 1} = wavelet_coeffs{lvl + 1}(unique_ctrl_pts_pointers_inds);
    reconstructed_control_points{lvl + 1} = reconstructed_control_points{lvl + 1}(unique_ctrl_pts_pointers_inds);
    
    if debug_flag == 1
        wcfs_time_lvl = toc(start_wcfs_time_lvl);
        disp(['Time taken to compute wavelet coefficients and reconstructed control points for level ' num2str(lvl + 1) ': ' num2str(wcfs_time_lvl) ' seconds']);
        disp('------------------------------------------------------------');
    end
    %profile viewer
end %End lvl
wcfs_time = toc(start_wcfs_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to compute all wavelet coefficients and reconstructed control points: ' num2str(wcfs_time) ' seconds']);
disp('************************************************************');

%profile off