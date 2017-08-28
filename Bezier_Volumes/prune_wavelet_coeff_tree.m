function pruned_wavelet_coeffs = prune_wavelet_coeff_tree(wavelet_coeffs, toprune, toprune2, shared_cells, ctrl_pts_pointers, myOT, start_lvl, max_lvl, b)

%Create a cell array that will contain the indices of the wavelet 
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
    %If we are not at the voxel level, we can look at the toprune2 cell
    %array, because toprune2 marks all octree cells that will be pruned
    %away (use this one instead of toprune, because toprune2 does not
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
                if isempty(shared_cells{lvl})||((~isempty(shared_cells{lvl}))&&(isempty(shared_cells{lvl, c})))
                    shared_cells{lvl, c} = ceil(find((c - ctrl_pts_pointers{lvl}) == 0)./8);
                end
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
            if isempty(shared_cells{lvl})||((~isempty(shared_cells{lvl}))&&(isempty(shared_cells{lvl, c})))
                shared_cells{lvl, c} = ceil(find((c - ctrl_pts_pointers{lvl}) == 0)./8);
            end
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
            if ~isempty(toprune{lvl - 1})
                test = find(ismember(parents, toprune{lvl - 1})) == 1;
                if length(test) == length(parents)
                    wc_toprune{lvl}(wc_toprune_cntr) = c;
                    wc_toprune_cntr = wc_toprune_cntr + 1;
                end  
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
