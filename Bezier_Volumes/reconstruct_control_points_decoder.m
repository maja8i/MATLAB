function reconstruction_decoder = reconstruct_control_points_decoder(debug_flag, rec_ctrlpts_forDec, wavelet_coeffs_forDec, SpatialIndex, FirstChildPtr, ChildCount, corner_coords_decoder, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, prune_flag, varargin)

%Assign inputs to new variables
rec_ctrlpts_start_lvl = rec_ctrlpts_forDec;
wavelet_coeffs = wavelet_coeffs_forDec;
if (~isempty(varargin)) && (prune_flag == 1)
    post_pruning_array = varargin{1};
    pp_first_nonempty = varargin{2};    %First octree level at which leaf cells are found (after pruning)
end

start_ctrlpts_recon_time = tic;

%Find the first octree level at which SpatialIndex is empty (assume that
%SpatialIndex is also empty at all levels higher than first_SI_empty)
first_SI_empty = find(cellfun(@isempty, SpatialIndex), 1);

%Initialize a cell array that will contain, for each octree level, the
%indices of the corners at this level that do not have wavelet
%coefficients associated with them (they are associated with low-pass
%coefficients instead)
cnrs_to_discard_all = cell((b + 1), 1);
%Initialize an averages cell array that will be used to store the average
%parent control points for each child corner for each occupied octree cell, 
%at each octree level (only for the corners that are associated with
%wavelet coefficients)
averages_all = cell((b + 1), 1);
%Initialize a cell array to store the indices of corners for each occupied
%octree cell at each octree level, which are associated with wavelet
%coefficients (these correspond to the "averages_all" cell array, above)
cnr_coords_inds_all = cell((b + 1), 1);
%Initialize a cell array that will contain the corresponding old indices of
%remaining corners after the corners that do not have wavelet coefficients
%associated with them have been removed, at each octree level
old_inds = cell((b + 1), 1);
%The below cell array will contain the ctrl_pts_pointers to unique corners
%that contain wavelet coefficients (and not low-pass coefficients)
ctrl_pts_pointers_wavelets = cell((b + 1), 1);
%The below cell array will contain unscaled ctrl_pts_pointers wavelets
unscaled_ctrl_pts_pointers_wavelets = ctrl_pts_pointers;

%Initialize a cell array to store the reconstructed signal (control point)
%at each vertex of each octree cell at every level from start_lvl to the 
%leaves
reconstruction_decoder = cell((b + 1), 1);

%Dequantize the reconstructed control points for start_lvl, which were 
%received from the encoder, then place them in the corresponding cell of 
%reconstruction_decoder. Note that when q_stepsize = 1, the dequantized
%values will just be equal to the quantized values.
reconstruction_decoder{start_lvl} = dequantize_uniform_scalar(rec_ctrlpts_start_lvl, q_stepsize);

%Start from the start_lvl and work up to the level before the voxel level
%(don't need to reconstruct the SDF at the voxel level, since it will not
%be used for anything there) ...
if isempty(first_SI_empty)
    %end_lvl = max_lvl - 1;
    end_lvl = max_lvl - 2;
else
    end_lvl = first_SI_empty - 2;
end

for lvl = start_lvl:end_lvl
    if debug_flag == 1
        disp(['Reconstructing control points at octree level ' num2str(lvl + 1) ' ...']);
        start_cp_recon_time_lvl = tic;
    end
%     %Dequantize the wavelet coefficients for all the corners (not just the
%     %unique ones) of all the occupied children (at level lvl + 1) of all 
%     %the occupied octree cells at the current level (lvl) - that is, get
%     %the dequantized wavelet coefficients for all the occupied octree cells
%     %at level lvl + 1, which have wavelet coefficients associated with them
%     all_wavelet_coeffs = dequantize_uniform_scalar(wavelet_coeffs{lvl + 1}(ctrl_pts_pointers{lvl + 1}), q_stepsize);
    %Counter for internal occupied octree cells at each level
    oc_code_cntr = 0;
    %For each occupied octree cell at the current level, if this cell has
    %children (i.e., is not a leaf), reconstruct the control points of all
    %the corners of all the cell's children
    for occ_cell = 1:size(SpatialIndex{lvl}, 1)
        %If this cell is a leaf
        if (prune_flag == 1) && ((lvl >= pp_first_nonempty) && (post_pruning_array{lvl}(occ_cell) == 1))
            %It has no children (they were pruned off at the encoder), so
            %do nothing further
            continue;  
        end
        %If this cell is not a leaf (i.e., it is internal), increment
        %oc_code_cntr
        oc_code_cntr = oc_code_cntr + 1;
        %Extract the cell's 8 corner coordinates. This cell will represent
        %our parent cell at the current level, since we will consider its
        %children for the control point (signal) reconstruction.
        parent_corner_coords = corner_coords_decoder{lvl}(((occ_cell*8 - 7):(occ_cell*8)), :); 
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
        child_ptrs = FirstChildPtr{lvl}(oc_code_cntr):(FirstChildPtr{lvl}(oc_code_cntr) + uint32(ChildCount{lvl}(oc_code_cntr)) - 1);
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
%         %In this case, the signal on each of these corners is a low-pass 
%         %coefficient (not a wavelet coefficient) and has already been 
%         %reconstructed as it is equal to its corresponding parent control 
%         %point. So, do nothing but place the value of the parent control 
%         %point in the corresponding location in the "averages" array - this 
%         %will ensure that the wavelet coefficient for the corresponding 
%         %child will be equal to 0 and the reconstructed control point for 
%         %the child will be the same as the parent control point.
%         if ~isempty(on_parent_child_inds)
%             averages(on_parent_child_inds) = parent_control_points(section_nbr);
%         end
        %In this case, the signal on each of these corners is a low-pass 
        %coefficient (not a wavelet coefficient) and has already been 
        %reconstructed, as it is equal to its corresponding parent control 
        %point
        if ~isempty(on_parent_child_inds)
            %reconstruction_decoder{lvl + 1}(on_parent_child_inds, 1) = parent_control_points(section_nbr);
            reconstruction_decoder{lvl + 1}(cnr_coords_inds(on_parent_child_inds), 1) = parent_control_points(section_nbr);
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

        %Remove the corners that do not have wavelet coefficients
        %associated with them (the control points for these corners have
        %already been reconstructed)
        cnrs_to_discard_all{lvl + 1}((end + 1):(end + length(on_parent_child_inds)), 1) = cnr_coords_inds(on_parent_child_inds);
        cnr_coords_inds(on_parent_child_inds) = [];
        averages(on_parent_child_inds) = [];

        cnr_coords_inds_all{lvl + 1}((end + 1):(end + length(cnr_coords_inds)), 1) = cnr_coords_inds;
        averages_all{lvl + 1}((end + 1):(end + length(cnr_coords_inds)), 1) = averages;
             
%         %Get the child cell wavelet coefficients for only the corners of 
%         %the child cells that have wavelet coefficients associated with them 
%         current_wavelet_coeffs = wavelet_coeffs{lvl + 1}(ctrl_pts_pointers{lvl + 1}(wav_cf_cntr:(wav_cf_cntr + length(cnr_coords_inds) - 1), 1));
% 
%         %Dequantize current_wavelet_coeffs
%         child_wavelet_coeffs = dequantize_uniform_scalar(current_wavelet_coeffs, q_stepsize);
%         
%         %Add the dequantized wavelet coefficients above, to the 
%         %corresponding values in "averages", to obtain the reconstructed 
%         %control point at each corresponding child corner
%         reconstruction_decoder{lvl + 1}(cnr_coords_inds, 1) = child_wavelet_coeffs + averages;    
         
    end %End occ_cell
    
    %For the corners that will remain at lvl + 1 after discarding those
    %corners in cnrs_to_discard_all, figure out what their corresponding
    %old indices would have been (before discarding the other corners)
    all_inds = 1:length(ctrl_pts_pointers{lvl + 1});
    old_inds{lvl + 1} = all_inds(ismember(all_inds, cnrs_to_discard_all{lvl + 1}) == 0);
    %Remove the pointers to the corners that do not have wavelet
    %coefficients associated with them
    unscaled_ctrl_pts_pointers_wavelets{lvl + 1}(cnrs_to_discard_all{lvl + 1}) = [];
    %Adjust the pointers to the wavelet_coeffs cell array, after corners
    %that are not associated with wavelet coefficients have been discarded,
    %i.e., scale the modified ctrl_pts_pointers_wavelets at the current 
    %octree level, to span only the remaining unique corners.
    [~, ~, unique_remaining_cnr_inds] = unique(unscaled_ctrl_pts_pointers_wavelets{lvl + 1}, 'stable');
    ctrl_pts_pointers_wavelets{lvl + 1} = unique_remaining_cnr_inds;

    %For each corner index in cnr_coords_inds_all at lvl + 1, get its new
    %index after the corners without wavelet coefficients have been
    %discarded
    new_inds = find(ismember(old_inds{lvl + 1}, cnr_coords_inds_all{lvl + 1}));

    %Get the (quantized) wavelet coefficients for ALL the corners of the 
    %child cells at lvl + 1 that have wavelet coefficients associated with 
    %them (not just the unique corners) 
    current_wavelet_coeffs = wavelet_coeffs{lvl + 1}(ctrl_pts_pointers_wavelets{lvl + 1}(new_inds));

    %Dequantize the quantized wavelet coefficients found above
    child_wavelet_coeffs = dequantize_uniform_scalar(current_wavelet_coeffs, q_stepsize);
        
    %Add the dequantized wavelet coefficients above, to their corresponding 
    %values in "averages", to obtain the reconstructed control point at 
    %each corresponding child corner that has a wavelet coefficient
    %associated with it
    reconstruction_decoder{lvl + 1}(cnr_coords_inds_all{lvl + 1}, 1) = child_wavelet_coeffs + averages_all{lvl + 1};
    
    %Keep only the reconstructed control points for the UNIQUE child 
    %corners at level lvl + 1, and discard the rest
    [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers{lvl + 1}, 'stable');
    reconstruction_decoder{lvl + 1} = reconstruction_decoder{lvl + 1}(unique_ctrl_pts_pointers_inds, 1);

    if debug_flag == 1
        cp_recon_time_lvl = toc(start_cp_recon_time_lvl);
        disp(' ');
        disp(['Time taken to reconstruct control points at level ' num2str(lvl + 1) ': ' num2str(cp_recon_time_lvl) ' seconds']);
    
        same_sign_cntr = 0;
        zero_cp_cntr = 0;
        for occ_cell = 1:size(SpatialIndex{lvl + 1}, 1)
            %Get all 8 control points for the corners of the current cell 
            current_ctrlpts = reconstruction_decoder{lvl + 1}(ctrl_pts_pointers{lvl + 1}((occ_cell*8 - 7):(occ_cell*8)));
            %Check if all control points of the current cell have the same
            %sign, including the case where all the control points may be 0
            if (abs(sum(sign(current_ctrlpts))) == 8)||(~any(sign(current_ctrlpts)))
                same_sign_cntr = same_sign_cntr + 1;
                %disp(['Cell ' num2str(occ_cell) ' has all control points with the same sign: ']);
                if ~any(sign(current_ctrlpts))
                    zero_cp_cntr = zero_cp_cntr + 1;
                end
            end
        end %End occ_cell
        disp(' ');
        disp(['TOTAL number of cells with all control points having the same sign, at level ' num2str(lvl + 1) ': ' num2str(same_sign_cntr)]);
        disp(['No. of cells with all 0 control points at level ' num2str(lvl + 1) ': ' num2str(zero_cp_cntr)]);
        disp('------------------------------------------------------------');
    end   
end %End lvl
ctrlpts_recon_time = toc(start_ctrlpts_recon_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to reconstruct all control points: ' num2str(ctrlpts_recon_time) ' seconds']);
disp('************************************************************');