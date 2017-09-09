function [pruned_occupancy_codes, post_pruning_array, toprune, toprune2] = prune_octree_old(myOT, all_zero_wav_cfs, start_lvl, max_lvl, b)

%Create a copy of the myOT.OccupancyCode cell array, which will contain
%the pruned set of occupancy codes
pruned_occupancy_codes = myOT.OccupancyCode;
%Create a cell array that will contain the indices of the octree cells 
%whose occupancy codes (which describe the locations of their occupied 
%children) are to be pruned away at each octree level, from
%pruned_occupancy_codes. These will be the internal cells only, not the
%leaves.
toprune = cell(size(myOT.OccupancyCode));
%Create a supplementary cell array that will contain a 1 where the 
%corresponding occupied cell is a leaf node (after pruning) and a 0 
%where the corresponding occupied cell is internal (not a leaf node).
post_pruning_array = cell(size(myOT.OccupancyCode));
%For each octree level that contains wavelet coefficients, except the 
%voxel level (since this level does not have its own occupancy codes,
%as the voxels do not have children), initialize all locations in 
%post_pruning_array to 0 (the leaf cells ("1" locations) will be
%determined after pruning)
for lvl = (start_lvl + 1):(max_lvl - 1)
    %If max_lvl is the voxel level (b + 1), then post_pruning_array
    %will only have (max_lvl - 1) levels, as these are the levels that
    %contain occupancy codes; if max_lvl is NOT b + 1, then the last
    %level in post_pruning_array will be set to 1 bits anyway (see
    %section below)
    post_pruning_array{lvl} = zeros(length(myOT.OccupancyCode{lvl}), 1);
end   
%If max_lvl is not the voxel level (b + 1), then set all the bits in 
%the last level in post_pruning_array to 1
if max_lvl < b + 1
     post_pruning_array{max_lvl} = ones(length(myOT.OccupancyCode{max_lvl}), 1);
end   
%Create a cell array that will contain the indices of the octree blocks
%that need to be pruned away from post_pruning_array, to account for 
%the internal nodes that were removed in pruned_occupancy_codes.
%NOTE: pruned_occupancy_codes and post_pruning_array will not be the 
%same size at the end: pruned_occupancy_codes will contain one 
%occupancy code per INTERNAL octree cell at each octree level, whereas 
%post_pruning_array will contain a 1 or 0 bit for EACH of the remaining 
%octree cells after pruning, including the internal cells AND the 
%leaves (except the leaves at level b + 1).
toprune2 = cell(size(myOT.OccupancyCode));
%Initialize a cell array that will keep track of the parent cells at lvl2
%levels (below) that have already been processed, to avoid repeating the
%same process for different children of the same parent
parents_already_processed = cell((b - 1), 1);

%The below code prunes branches (occupancy codes) of octree cells that 
%have all zero wavelet coefficients throughout the branch, so the 
%pruning is done bottom-up from the voxel level ...
tic;
lvl = max_lvl; %Here assume max_lvl = b + 1
%If there exist any occupied voxels, which contain all zero wavelet 
%coefficients ...
if ~isempty(all_zero_wav_cfs{lvl})
    %disp('Found voxels with all-zero wavelet coefficients: ');
    %disp(num2str(all_zero_wav_cfs{lvl}'));
    %Get the parent cells for all of these voxels
    parents = myOT.ParentPtr{lvl}(all_zero_wav_cfs{lvl});
    %Extract only the unique parent cells, since the parent indices in
    %"parents", above, will be repeated for all the children of a parent
    parents = unique(parents, 'stable');
    %Find the first child pointers for each of the parents, and convert
    %their child counts into uint32 type (same as the child pointers), so
    %that these can be added to the child pointers to get the last child
    %pointers
    first_child = myOT.FirstChildPtr{lvl - 1}(parents);
    child_count = uint32(myOT.ChildCount{lvl - 1}(parents));
    last_child = first_child + child_count - 1;
    %For each of the unique parents of the all-zero-wavelet-coefficient
    %voxels ...
    for p_cell = 1:length(parents)
        %profile on
        %Get the current parent
        parent_cell = parents(p_cell);
        %Check if this parent_cell has already been marked for pruning (it
        %may have been a child of a parent whose occupancy code was marked
        %for pruning in a previous iteration)
        if any(toprune{lvl - 1} == parent_cell)
            %No need to check this parent_cell again
            continue;
        end
        %Find all the children (voxels) of parent_cell 
        children = first_child(p_cell):last_child(p_cell);
        %Check if ALL the children (voxels) have all zero wavelet
        %coefficients
        az_children = zeros(1, length(children));
        for i = 1:length(children)
            az_children(i) = any(all_zero_wav_cfs{lvl} == children(i));
        end
        %If ALL of the child voxels have all zero wavelet coefficients
        if all(az_children)
            %Mark the occupancy code of parent_cell for pruning
            toprune{lvl - 1}((length(toprune{lvl - 1}) + 1), 1) = parent_cell;
            %If the current parent_cell does NOT have all zero wavelet 
            %coefficients (but all of its child voxels do)              
            if ~any(all_zero_wav_cfs{lvl - 1} == parent_cell)
                %The current parent_cell will become a leaf after pruning 
                %its children away
                post_pruning_array{lvl - 1}(parent_cell) = 1;
                %Do not need to check ancestors of the current parent_cell, 
                %so move on to checking the next parent_cell instead
                %disp(['Finished processing parent_cell ' num2str(p_cell) '/' num2str(length(parents))]);
                continue;
            end
            %If the parent_cell does have all zero wavelet coefficients (as
            %do all of its children), successively check ancestors of the
            %parent_cell, to see if we can prune further up the octree 
            for lvl2 = (lvl - 1):-1:(start_lvl + 1)
                %Get the parent of the current parent_cell
                parent2 = myOT.ParentPtr{lvl2}(parent_cell);
                %Check if parent2 has already been processed in another
                %iteration (for a different child (parent_cell), perhaps)
                %or if it has already been marked for pruning in a
                %different iteration
                if (~isempty(parents_already_processed{lvl2 - 1}) && any(parents_already_processed{lvl2 - 1} == parent2)) || (any(toprune{lvl2 - 1} == parent2))
                    break;
                end
                %Record parent2 as being processed, so that it will not be
                %repeated in a future iteration
                parents_already_processed{lvl2 - 1}(length(parents_already_processed{lvl2 - 1}) + 1) = parent2;
                %Get the children of parent2 - one of these children will 
                %be the current parent_cell
                children2 = ((myOT.FirstChildPtr{lvl2 - 1}(parent2)):(myOT.FirstChildPtr{lvl2 - 1}(parent2) + uint32((myOT.ChildCount{lvl2 - 1}(parent2))) - 1));
                %Check if children2 have all zero wavelet coefficients
                az_children2 = zeros(1, length(children2));
                for i = 1:length(children2)
                    az_children2(i) = any(all_zero_wav_cfs{lvl2} == children2(i));
                end
                %If NOT all of parent2's children (parent_cell's siblings) 
                %have all zero wavelet coefficients
                if ~all(az_children2)
                    %We cannot prune the occupancy code of parent2, but the 
                    %current parent_cell will become a leaf
                    post_pruning_array{lvl2}(parent_cell) = 1;
                    %Stop searching other ancestors of the current 
                    %parent_cell; move on to the next parent_cell
                    break;
                end
                %If ALL of parent2's children have all zero wavelet
                %coefficients, then the current parent_cell will not be a 
                %leaf, so mark it for pruning from post_pruning array
                toprune2{lvl2}((length(toprune2{lvl2}) + 1), 1) = parent_cell;
                %Mark the occupancy code for the current ancestor parent 
                %for pruning
                toprune{lvl2 - 1}((length(toprune{lvl2 - 1}) + 1), 1) = parent2;
                %Mark the occupancy codes of parent2's other occupied 
                %children (i.e., parent_cell's siblings) for pruning too 
                %(since these children are occupied, they must have 
                %occupancy codes for children of their own). One of the 
                %children2 will be the current parent_cell, which has
                %already been included in toprune and toprune2, so we 
                %exclude it in the lines below.
                toprune{lvl2}((length(toprune{lvl2}) + 1):(length(toprune{lvl2}) + length(children2) - 1), 1) = children2(children2 ~= parent_cell);
                %These children will not be leaves (since they will be 
                %pruned off), so mark them for pruning from 
                %post_pruning_array 
                toprune2{lvl2}((length(toprune2{lvl2}) + 1):(length(toprune2{lvl2}) + length(children2) - 1), 1) = children2(children2 ~= parent_cell);
                %If parent2 does NOT have all zero wavelet coefficients 
                %(but all of its children do)
                if ~any(all_zero_wav_cfs{lvl2 - 1} == parent2)
                    %parent2 will become a leaf
                    post_pruning_array{lvl2 - 1}(parent2) = 1; 
                    %Stop checking other ancestors
                    break;
                end
                %If parent2 does have all zero wavelet coefficients (as do 
                %all of its children), continue checking other ancestors, 
                %with the current parent2 being the new parent_cell
                parent_cell = parent2; 
            end %End for-loop for lvl2
        %If NOT ALL of the children (voxels) of the current parent_cell
        %have all zero wavelet coefficients, we cannot prune the occupancy 
        %code of the current parent_cell, so move on to checking the next 
        %parent_cell
        end %End check if all(az_children)     
        %profile viewer
        %disp(['Finished processing parent_cell ' num2str(p_cell) '/' num2str(length(parents))]);
    end %End p_cell
end %End check if ~isempty(all_zero_wav_cfs{lvl})

%For each octree level covered in toprune ...
for lvl = 1:size(toprune, 1)
    if ~isempty(toprune{lvl})
        %Prune the occupancy codes of the octree cells that have been 
        %marked
        pruned_occupancy_codes{lvl}(toprune{lvl}) = [];
        disp(['Level ' num2str(lvl) ' octree pruning done (removed occupancy codes of leaf cells at this level, and their descendants)']);
        disp('------------------------------------------------------------');
    end
end

%For each octree level covered in toprune2 ...
for lvl = 1:size(toprune2, 1)
    if ~isempty(toprune2{lvl})
        %Prune the post_pruning_array accordingly
        post_pruning_array{lvl}(toprune2{lvl}) = [];
        disp(['Level ' num2str(lvl) ' post-pruning array pruning done']);
        disp('------------------------------------------------------------');
    end
end

ot_pruning_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to prune octree: ' num2str(ot_pruning_time) ' seconds']);
disp('************************************************************');