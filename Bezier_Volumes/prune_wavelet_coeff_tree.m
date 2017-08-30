function pruned_wavelet_coeffs = prune_wavelet_coeff_tree(wavelet_coeffs, toprune, toprune2, ctrl_pts_pointers, myOT, start_lvl, max_lvl, b)

%Create a copy of the wavelet_coeffs cell array, which will contain the 
%pruned set of wavelet coefficients
pruned_wavelet_coeffs = wavelet_coeffs;

%For each octree level at which there are wavelet coefficients ...
tic;
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
            %Remove the wavelet coefficients associated with the corners of
            %the pruned_cells at this octree level
            corners_to_prune = ((pruned_cells(1)*8) - 7):(pruned_cells(1)*8);
            for i = 2:length(pruned_cells)
                corners_to_prune = cat(1, corners_to_prune, ((pruned_cells(i)*8) - 7):(pruned_cells(i)*8));
            end
            wavelet_coeffs_expanded(corners_to_prune) = [];
            %Remove the corresponding corners in ctrl_pts_pointers_temp
            ctrl_pts_pointers_temp(corners_to_prune) = [];
            %Keep only the wavelet coefficients corresponding to the unique
            %corners that remain after pruning
            pruned_wavelet_coeffs{lvl} = wavelet_coeffs_expanded(unique(ctrl_pts_pointers_temp, 'stable'));
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
            %Find the children (voxels) of pruned_parent_cells
            children_cnt = 1;
            for j = 1:length(pruned_parent_cells)
                children_of_pruned(children_cnt:(children_cnt + myOT.ChildCount{lvl - 1}(pruned_parent_cells(j)) - 1)) = (myOT.FirstChildPtr{lvl - 1}(pruned_parent_cells(j))):(myOT.FirstChildPtr{lvl - 1}(pruned_parent_cells(j)) + uint32(myOT.ChildCount{lvl - 1}(pruned_parent_cells(j))) - 1);
                children_cnt = children_cnt + myOT.ChildCount{lvl - 1}(pruned_parent_cells(j));
            end
            %Expand wavelet_coeffs{lvl} to get ALL the corners at this
            %octree level, not just the unique ones
            wavelet_coeffs_expanded = wavelet_coeffs{lvl}(ctrl_pts_pointers{lvl});
            %Make a copy of ctrl_pts_pointers{lvl}, as it will be
            %temporarily modified
            ctrl_pts_pointers_temp = ctrl_pts_pointers{lvl};
            %Remove the wavelet coefficients associated with the corners of
            %the children_of_pruned voxels at this octree level
            corners_to_prune = ((children_of_pruned(1)*8) - 7):(children_of_pruned(1)*8);
            for i = 2:length(children_of_pruned)
                corners_to_prune = cat(1, corners_to_prune, ((children_of_pruned(i)*8) - 7):(children_of_pruned(i)*8));
            end
            wavelet_coeffs_expanded(corners_to_prune) = [];
            %Remove the corresponding corners in ctrl_pts_pointers_temp
            ctrl_pts_pointers_temp(corners_to_prune) = [];
            %Keep only the wavelet coefficients corresponding to the unique
            %corners that remain after pruning
            pruned_wavelet_coeffs{lvl} = wavelet_coeffs_expanded(unique(ctrl_pts_pointers_temp, 'stable'));
        end %End check if ~isempty(toprune{lvl - 1})   
    end %End check if lvl < (b + 1)
    disp(['Level ' num2str(lvl) ' done']);
    disp('------------------------------------------------------------');
    %profile viewer
end %End lvl   
%profile off

w_pruning_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to prune wavelet tree: ' num2str(w_pruning_time) ' seconds']);
disp('************************************************************');