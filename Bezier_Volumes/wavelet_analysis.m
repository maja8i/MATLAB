function [wavelet_coeffs, reconstructed_control_points, all_zero_wav_cfs, ctrl_pts_pointers_wavelets, cnrs_to_discard_all, old_inds] = wavelet_analysis(debug_flag, myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, zero_threshold)

%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) that have 
%wavelet coefficients associated with them (and not low-pass coefficients),
%across all octree blocks and levels, starting from start_lvl and going up 
%to one level before the leaves
wavelet_coeffs = cell((b + 1), 1);  %These will be quantized coefficients
%Initialize a cell array to store the reconstructed signal at each vertex
%of each octree cell at every level from start_lvl to one level before the
%leaves
reconstructed_control_points = cell(size(control_points, 1), 1);
%The below cell array will contain the ctrl_pts_pointers to unique corners
%that contain wavelet coefficients (and not low-pass coefficients)
ctrl_pts_pointers_wavelets = cell((b + 1), 1);
%The below cell array will contain unscaled ctrl_pts_pointers wavelets
unscaled_ctrl_pts_pointers_wavelets = ctrl_pts_pointers;
%Initialize a cell array to store the indices of the occupied octree cells
%at each level that have quantized wavelet coefficient values (i.e., 
%quantized symbols) of 0 on all of their corners that are associated with
%wavelet coefficients
all_zero_wav_cfs = cell((b + 1), 1);
%Initialize a cell array that will contain, for each octree level, the
%indices of the corners at this level that do not have wavelet
%coefficients associated with them (they are associated with low-pass
%coefficients instead)
cnrs_to_discard_all = cell((b + 1), 1);
%Initialize a cell array that will contain the corresponding old indices of
%remaining corners after the corners that do not have wavelet coefficients
%associated with them have been removed, at each octree level
old_inds = cell((b + 1), 1);

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
    %Initialize a counter for the wavelet coefficients at the current
    %octree level
    wav_cf_cntr = 1;
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        if (lvl == 5)&&(occ_cell == 174)
            pause(1);
        end
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
%         %In this case, the signal on each of these corners is a low-pass 
%         %coefficient (not a wavelet coefficient) and has already been 
%         %reconstructed as it is equal to its corresponding parent control 
%         %point. So, do nothing but place the value of the reconstructed 
%         %parent control point in the corresponding location in the
%         %"averages" array and change the value of the corresponding child 
%         %control point to be the same as the value of the reconstructed 
%         %parent control point - this will ensure that the wavelet 
%         %coefficient for the corresponding child will be equal to 0 and the
%         %reconstructed control point for the child will be the same as the
%         %parent control point.
%         if ~isempty(on_parent_child_inds)
%             averages(on_parent_child_inds) = parent_control_points(section_nbr);
%             child_control_points(on_parent_child_inds) = parent_control_points(section_nbr);
%         end
        %In this case, the signal on each of these corners is a low-pass 
        %coefficient (not a wavelet coefficient) and has already been 
        %reconstructed, as it is equal to its corresponding parent control 
        %point
        if ~isempty(on_parent_child_inds)
            %reconstructed_control_points{lvl + 1}(on_parent_child_inds, 1) = parent_control_points(section_nbr);
            reconstructed_control_points{lvl + 1}(cnr_coords_inds(on_parent_child_inds), 1) = parent_control_points(section_nbr);
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
        
        %For all the child corners of the current parent block, except 
        %those that lie on parent corners, subtract their corresponding 
        %average parent signal (control point) computed above, from the 
        %child control point. The result will be the high-pass transform 
        %(wavelet) coefficient of the child corner.
        cnrs_to_discard_all{lvl + 1}((end + 1):(end + length(on_parent_child_inds)), 1) = cnr_coords_inds(on_parent_child_inds);
        cnr_coords_inds(on_parent_child_inds) = [];
        child_control_points(on_parent_child_inds) = [];
        averages(on_parent_child_inds) = [];
        %wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = child_control_points - averages;
        wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1) = child_control_points - averages;     
        %Quantize the wavelet coefficients computed above
        %wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(cnr_coords_inds, 1), q_stepsize);
        wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1), q_stepsize);
        %Check if any of the quantized wavelet coefficients are within +/-
        %zero_threshold of 0: if so, set the values of these quantized
        %coefficients to 0
        thresholded_wav_cfs = wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1);
        thresholded_wav_cfs(abs(thresholded_wav_cfs) <= zero_threshold) = 0;
        wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1) = thresholded_wav_cfs;
        %Check if all the wavelet coefficients for any of the child cells
        %of the current occ_cell are zero (or within zero_threshold), and 
        %if so then record these child cells' indices, to consider them for
        %pruning later
        current_wavelet_coeffs = wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1);
        cells_with_parent_corners = ceil(on_parent_child_inds./8);
        wcf_cntr = 1;
        for c_cell = 1:length(child_ptrs)
            shared_parent_cnrs = find(cells_with_parent_corners == c_cell);
            %If any of the current child cell corners are shared with a
            %parent corner
            if ~isempty(shared_parent_cnrs)
                %Only get the wavelet coefficients for the other
                %corners of this child cell
                child_wav_cfs = current_wavelet_coeffs(wcf_cntr:(wcf_cntr + 7 - length(shared_parent_cnrs)));
                wcf_cntr = wcf_cntr + 8 - length(shared_parent_cnrs);
            else
                %Get the wavelet coefficients for all 8 corners of this
                %child cell
                child_wav_cfs = current_wavelet_coeffs(wcf_cntr:(wcf_cntr + 7));
                wcf_cntr = wcf_cntr + 8;     
            end
%                 %If the quantized wavelet coefficients at all the corners of this 
%                 %cell are 0 ...
%                 if isempty(find((current_wavelet_coeffs ~= 0), 1))
            %If all the quantized wavelet coefficients for the current
            %child cell are near 0 (i.e., within +/- zero_threshold of
            %0) ...
            if isempty(find((abs(child_wav_cfs) > zero_threshold), 1))
                all_zero_wav_cfs{lvl + 1}(end + 1) = child_ptrs(c_cell); 
            end
        end
        
        %Dequantize the quantized wavelet coefficients and add them to the
        %corresponding values in "averages", to obtain the reconstructed 
        %signal (control point) at each corresponding child corner. Note
        %that the QUANTIZED versions of these wavelet coefficients will be
        %sent to the decoder.
        dequant_wavelet_coeffs = dequantize_uniform_scalar(wavelet_coeffs{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1), q_stepsize);
        reconstructed_control_points{lvl + 1}(cnr_coords_inds, 1) = averages + dequant_wavelet_coeffs;
        
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
        
        wav_cf_cntr = wav_cf_cntr + length(cnr_coords_inds);
        
        %disp(['Finished occ_cell ' num2str(occ_cell) '/' num2str(myOT.NodeCount(lvl))]);
    end %End occ_cell
    
    %Keep only the wavelet coefficients and reconstructed control points
    %for the UNIQUE child corners at level lvl + 1, and discard the rest
    [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers{lvl + 1}, 'stable');
    reconstructed_control_points{lvl + 1} = reconstructed_control_points{lvl + 1}(unique_ctrl_pts_pointers_inds);
    %For the corners that will remain at lvl + 1 after discarding those
    %corners in cnrs_to_discard_all, figure out what their corresponding
    %old indices would have been (before discarding the other corners)
    all_inds = 1:length(unscaled_ctrl_pts_pointers_wavelets{lvl + 1});
    old_inds{lvl + 1} = all_inds(ismember(all_inds, cnrs_to_discard_all{lvl + 1}) == 0);
    unscaled_ctrl_pts_pointers_wavelets{lvl + 1}(cnrs_to_discard_all{lvl + 1}) = [];  %Remove the discarded corners' (those that don't have wavelet coefficients) pointers
    [~, unique_ctrl_pts_pointers_wavelets_inds, unique_remaining_cnr_inds] = unique(unscaled_ctrl_pts_pointers_wavelets{lvl + 1}, 'stable');
    wavelet_coeffs{lvl + 1} = wavelet_coeffs{lvl + 1}(unique_ctrl_pts_pointers_wavelets_inds);
    %Scale the modified ctrl_pts_pointers_wavelets at the current octree
    %level, to span only the remaining unique corners
    ctrl_pts_pointers_wavelets{lvl + 1} = unique_remaining_cnr_inds;

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