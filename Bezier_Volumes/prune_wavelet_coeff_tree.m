function pruned_wavelet_coeffs = prune_wavelet_coeff_tree(debug_flag, wavelet_coeffs, toprune, toprune2, ctrl_pts_pointers, myOT, start_lvl, max_lvl, b)

%Create a copy of the wavelet_coeffs cell array, which will contain the 
%pruned set of wavelet coefficients
pruned_wavelet_coeffs = wavelet_coeffs;

%For each octree level at which there are wavelet coefficients ...
start_w_pruning_time = tic;
for lvl = (start_lvl + 1):max_lvl
    %profile on
    %If we are not at the voxel level, we can look at the toprune2 cell
    %array. Use toprune2 instead of toprune, because toprune2 does not mark
    %the leaf cells for pruning and we want to keep the wavelet 
    %coefficients for the leaf cells (since their corner coordinates will
    %be reconstructed at the decoder and we need to have a wavelet 
    %coefficient per reconstructed corner).
    if lvl < (b + 1)
        %If this level contains any octree cells whose occupancy codes have 
        %been pruned ...
        if ~isempty(toprune2{lvl})
            %Get the indices of these cells
            pruned_cells = toprune2{lvl};
            %Expand wavelet_coeffs{lvl} to get ALL the corners at this
            %octree level, not just the unique ones
            wavelet_coeffs_expanded = wavelet_coeffs{lvl}(ctrl_pts_pointers{lvl});
            %Make a copy of ctrl_pts_pointers{lvl}, as it will be
            %temporarily modified
            ctrl_pts_pointers_temp = ctrl_pts_pointers{lvl};
            %Find the corners to prune at this octree level 
            corners_to_prune = [(pruned_cells.*8 - 7) (pruned_cells.*8 - 6) (pruned_cells.*8 - 5) (pruned_cells.*8 - 4) (pruned_cells.*8 - 3) (pruned_cells.*8 - 2) (pruned_cells.*8 - 1) (pruned_cells.*8)];
            %Remove the wavelet coefficients associated with the 
            %corners_to_prune at this octree level
            wavelet_coeffs_expanded(corners_to_prune) = [];
            %Remove the corresponding corners in ctrl_pts_pointers_temp
            ctrl_pts_pointers_temp(corners_to_prune) = [];
            %Keep only the wavelet coefficients corresponding to the unique
            %corners that remain after pruning
            [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers_temp, 'stable');
            pruned_wavelet_coeffs{lvl} = wavelet_coeffs_expanded(unique_ctrl_pts_pointers_inds);
        end %End check if ~isempty(toprune2{lvl})
    %If lvl == (b + 1)
    else  
        %The voxel level is not covered in toprune2 or to prune, since 
        %voxels do not have occupancy codes of their own (as voxels have no 
        %children), but the voxels still have wavelet coefficients on their 
        %corners and some of these may be candidates for pruning. The 
        %voxel corners that can be pruned are only those that are 
        %descendants of octree cells that were pruned off at previous 
        %octree levels and are not shared with any non-pruned voxels. Look
        %inside toprune at the level just before the voxel level, to see if
        %any of the voxels' parents there were marked to have their
        %occupancy codes (and therefore children) pruned off ...
        if ~isempty(toprune{lvl - 1})
            pruned_parent_cells = toprune{lvl - 1};
            %Find the children (voxels) of each of the pruned_parent_cells
            first_child = myOT.FirstChildPtr{lvl - 1}(pruned_parent_cells);
            child_count = uint32(myOT.ChildCount{lvl - 1}(pruned_parent_cells));
            last_child = first_child + child_count - 1;
            children_cnt = 1;
            for j = 1:length(pruned_parent_cells)
                children_of_pruned(children_cnt:(children_cnt + child_count(j) - 1)) = first_child(j):last_child(j);
                children_cnt = children_cnt + child_count(j);
            end
            %Expand wavelet_coeffs{lvl} to get ALL the corners at this
            %octree level, not just the unique ones
            wavelet_coeffs_expanded = wavelet_coeffs{lvl}(ctrl_pts_pointers{lvl});
            %Make a copy of ctrl_pts_pointers{lvl}, as it will be
            %temporarily modified
            ctrl_pts_pointers_temp = ctrl_pts_pointers{lvl};
            %Find the corners to prune at this octree level
            corners_to_prune = [(children_of_pruned.*8 - 7) (children_of_pruned.*8 - 6) (children_of_pruned.*8 - 5) (children_of_pruned.*8 - 4) (children_of_pruned.*8 - 3) (children_of_pruned.*8 - 2) (children_of_pruned.*8 - 1) (children_of_pruned.*8)];
            %Remove the wavelet coefficients associated with the 
            %corners_to_prune at this octree level
            wavelet_coeffs_expanded(corners_to_prune) = [];
            %Remove the corresponding corners in ctrl_pts_pointers_temp
            ctrl_pts_pointers_temp(corners_to_prune) = [];
            %Keep only the wavelet coefficients corresponding to the unique
            %corners that remain after pruning
            [~, unique_ctrl_pts_pointers_inds, ~] = unique(ctrl_pts_pointers_temp, 'stable');
            pruned_wavelet_coeffs{lvl} = wavelet_coeffs_expanded(unique_ctrl_pts_pointers_inds);
        end %End check if ~isempty(toprune{lvl - 1})   
    end %End check if lvl < (b + 1)
    if debug_flag == 1
        disp(['Level ' num2str(lvl) ' done']);
        disp('------------------------------------------------------------');
    end
    %profile viewer
end %End lvl   
%profile off

w_pruning_time = toc(start_w_pruning_time);
disp(' ');
disp('************************************************************');
disp(['Time taken to prune wavelet coefficient tree: ' num2str(w_pruning_time) ' seconds']);
disp('************************************************************');
