function [pruned_occupancy_codes, post_pruning_array, toprune, toprune2] = prune_octree(myOT, all_zero_wav_cfs, start_lvl, max_lvl, b)

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
%Initialize a cell array that will keep track of the parent cells that have 
%already been processed at each octree level, to avoid repeating the same
%process for different children of the same parent
parents_already_processed = cell((b - 1), 1);

%The below code prunes branches (occupancy codes) of octree cells that 
%have all zero wavelet coefficients throughout the branch, so the 
%pruning is done bottom-up from the voxel level ...
tic;
lvl = max_lvl;  %Here assume max_lvl = b + 1
%If there exist any occupied voxels, which contain all zero wavelet 
%coefficients ...
if ~isempty(all_zero_wav_cfs{lvl})
    %Get the parent cells for all of these voxels
    parents = myOT.ParentPtr{lvl}(all_zero_wav_cfs{lvl});
    %Extract only the unique parent cells, since the parent indices 
    %in "parents", above, will be repeated for all the children of 
    %a parent
    parents = unique(parents, 'stable');
    %Find the first child pointers for each of the parents, and 
    %convert their child counts into uint32 type (same as the child
    %pointers), so that these can be added to the child pointers to 
    %get the last child pointers
    first_child = myOT.FirstChildPtr{lvl - 1}(parents);
    child_count = uint32(myOT.ChildCount{lvl - 1}(parents));
    last_child = first_child + child_count - 1;   
    %For each of the unique parents of the all-zero-wavelet-coeffi-
    %cient voxels ...
    for p_cell = 1:length(parents)
        %profile on
        %Get the current parent
        parent_cell = parents(p_cell);
        %Find all the children (voxels) of parent_cell 
        children = first_child(p_cell):last_child(p_cell);
        %Check if ALL the children (voxels) have all zero wavelet
        %coefficients
        az_children = zeros(1, length(children));
        for i = 1:length(children)
            az_children(i) = any(all_zero_wav_cfs{lvl} == children(i));
        end
        %If ALL of the child voxels have all zero wavelet 
        %coefficients
        if all(az_children)
            %Mark the occupancy code of parent_cell for pruning
            %(regardless of whether or not the parent_cell has all
            %zero wavelet coefficients itself)
            toprune{lvl - 1}((end + 1), 1) = parent_cell;
            %If the current parent_cell does NOT have all zero 
            %wavelet coefficients (but all of its child voxels do)              
            if ~any(all_zero_wav_cfs{lvl - 1} == parent_cell)
                %The current parent_cell will become a leaf after 
                %pruning its children away
                post_pruning_array{lvl - 1}(parent_cell) = 1;
            end
        end %End check if all(az_children)
        %disp(['Finished p_cell ' num2str(p_cell) '/' num2str(length(parents))]);
        %profile viewer
    end %End p_cell
end %End check if ~isempty(all_zero_wav_cfs{lvl})
%disp('Finished voxel level');

%Go through other octree levels, to see if we can prune further up the
%octree
for lvl = (max_lvl - 1):-1:(start_lvl + 1)   
    if ~isempty(toprune{lvl})
        %Check all of the octree cells marked in toprune for lvl
        for ot_cell = (toprune{lvl})' %Need to have row vector, not column vector, to use as indices into for-loop
            %If the cell was marked as a leaf, do not process it further 
            %(no need to check other ancestors); move on to the next 
            %ot_cell
            if post_pruning_array{lvl}(ot_cell) == 1
                continue;
            end
            %If the cell was not marked as a leaf, get its parent at 
            %lvl - 1
            parent = myOT.ParentPtr{lvl}(ot_cell);
            %Multiple cells at lvl may have the same parent, so check if 
            %parent has already been processed (we don't want to process it 
            %more than once)
            if (~isempty(parents_already_processed{lvl - 1}) && any(parents_already_processed{lvl - 1} == parent))
                continue;
            end
            %Record parent as being processed, so that it will not be 
            %repeated in a future iteration
            parents_already_processed{lvl - 1}(end + 1) = parent;
            %Get all the children of parent - one of these children will be 
            %the current ot_cell
            children = ((myOT.FirstChildPtr{lvl - 1}(parent)):(myOT.FirstChildPtr{lvl - 1}(parent) + uint32((myOT.ChildCount{lvl - 1}(parent))) - 1));
            %Check if all of these children have been marked for pruning 
            %(in the previous iteration)
            children_toprune = zeros(1, length(children));
            for i = 1:length(children)
                children_toprune(i) = any(toprune{lvl} == children(i));
            end
            %If NOT all of the children have been marked for pruning
            if ~all(children_toprune)
                %We cannot prune the occupancy code of parent, so the
                %current ot_cell and any of its siblings that have been
                %marked for pruning will become leaves
                post_pruning_array{lvl}(children(children_toprune == 1)) = 1;
                %Move on to checking the parent of the next ot_cell
                continue;
            end
            %If all of the children's occupancy codes have been marked for
            %pruning, check if any of the children were marked as leaf 
            %cells
            if any(post_pruning_array{lvl}(children))
                %If any of the children were marked as leaf cells, then we
                %cannot prune the occupancy code of the parent, so we need
                %to make all the children leaves
                post_pruning_array{lvl}(children) = 1;
                %Move on to checking the parent of the next ot_cell
                continue;
            end
            %If all of the children's occupancy codes were marked for 
            %pruning AND none of the children were marked as leaves, then
            %we can mark the current parent's occupancy code for pruning
            toprune{lvl - 1}((end + 1), 1) = parent;
            %Also mark the children for pruning from post_pruning_array, 
            %since none of them will be leaves (this will include the
            %current ot_cell, too)
            toprune2{lvl}((end + 1):(end + length(children)), 1) = children;
            %Check if the current parent has all zero wavelet coefficients: 
            %if not, but all of its children do (since they were all marked 
            %for pruning and none of them are leaves), then the parent will
            %be a leaf
            if ~any(all_zero_wav_cfs{lvl - 1} == parent)
                post_pruning_array{lvl - 1}(parent) = 1; 
            end      
        end %End ot_cell
        %disp(['Finished level ' num2str(lvl)]);
    end %End check if ~isempty(toprune{lvl})
end %End lvl

%For each octree level covered in toprune ...
for lvl = 1:size(toprune, 1)
    if ~isempty(toprune{lvl})
        %Prune the occupancy codes of the octree cells that have been 
        %marked
        pruned_occupancy_codes{lvl}(toprune{lvl}) = [];
        disp(['Level ' num2str(lvl) ' octree pruning done']);
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


