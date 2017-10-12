function [pruned_occupancy_codes, post_pruning_array, toprune, toprune2] = prune_octree(myOT, all_zero_wav_cfs, start_lvl, max_lvl, b, reconstructed_control_points, ctrl_pts_pointers)

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

%The below code prunes branches (occupancy codes) of octree cells that 
%have all zero wavelet coefficients throughout the branch, so the 
%pruning is done bottom-up from the voxel level ...
tic;
lvl = max_lvl;  %Here assume max_lvl = b + 1
%If there exist any occupied voxels, which contain all zero wavelet 
%coefficients ...
if ~isempty(all_zero_wav_cfs{lvl})
    %Get the parent cells for all of these voxels
    parents_all = myOT.ParentPtr{lvl}(all_zero_wav_cfs{lvl});
    %Extract only the unique parent cells, since the parent indices 
    %in parents_all, above, will be repeated for all the children of 
    %a parent
    [parents, ~, parents_pointers] = unique(parents_all, 'stable');
    %Get the child count for each of the unique parents and convert it to
    %int8 type (from uint8), in case the subtraction from unique_count,
    %below, gives a negative result (otherwise the negative result would
    %just be truncated to 0)
    child_count = int8(myOT.ChildCount{lvl - 1}(parents));
    %Find the number of times that each unique parent is referenced in
    %parents_all
    M = size(parents, 1);
    N = size(parents_all, 1); 
    A = sparse(parents_pointers, [1:N]', ones(N,1), M, N);
    unique_count = int8(full(sum(A, 2)));  %Convert to int8 to be the same type as child_count, so we can do the subtraction below
    %Find which of the unique parents has all of their children with all
    %zero wavelet coefficients: if the number of times that a unique parent
    %appears in parents_all is equal to the total number of children that 
    %this parent has, then all of this parent's children (voxels) have all 
    %zero wavelet coefficients. We only need to consider the parents of
    %these voxels further, since their occupancy codes can be pruned;
    %remove the other parents from the list, since we cannot prune their
    %occupancy codes.
    parents((unique_count - child_count) ~= 0) = [];
    %Mark the occupancy codes of all the remaining parents for pruning
    toprune{lvl - 1}((1:length(parents)), 1) = parents;
    %For each of the unique parents whose occupancy code has been marked
    %for pruning (i.e., those whose children all have all zero wavelet 
    %coefficients), check whether or not this parent will be a leaf after 
    %its children have been pruned away. At this stage, a parent will only 
    %be marked as a leaf if it does NOT have all zero wavelet coefficients
    %(but all of its children do, hence why they can be pruned away). 
    post_pruning_array{lvl - 1}(parents(ismember(parents, all_zero_wav_cfs{lvl - 1}) == 0)) = 1;
end %End check if ~isempty(all_zero_wav_cfs{lvl})
disp(['Finished marking level ' num2str(lvl) ' for pruning']);

%Go through other octree levels, to see if we can prune further up the
%octree
for lvl = (max_lvl - 1):-1:(start_lvl + 1)   
    %Only check the ancestors of octree cells that have been marked for
    %pruning at the current octree level
    if ~isempty(toprune{lvl})
        %Get the parents of the octree cells that have been marked for
        %pruning at the current octree level
        parents_all = myOT.ParentPtr{lvl}(toprune{lvl});
        %Extract only the unique parent cells from parents_all, since the 
        %parent indices in parents_all will be repeated for all of the
        %parents' children that are found in toprune{lvl}
        [parents, ~, parents_pointers] = unique(parents_all, 'stable');
        %Get the child count for each of the unique parents and convert it
        %to int8 type (from uint8), in case the subtraction from 
        %unique_count, below, gives a negative result (otherwise the 
        %negative result would just be truncated to 0)
        child_count = int8(myOT.ChildCount{lvl - 1}(parents));
        %Find the number of times that each unique parent is referenced in
        %parents_all
        M = size(parents, 1);
        N = size(parents_all, 1); 
        A = sparse(parents_pointers, [1:N]', ones(N,1), M, N);
        unique_count = int8(full(sum(A, 2)));  %Convert to int8 to be the same type as child_count, so we can do the subtraction below
        %Find which of the unique parents have all of their children's
        %occupancy codes marked for pruning in toprune{lvl}: if the number 
        %of times that a unique parent appears in parents_all is equal to 
        %the total number of children that this parent has, this means that
        %all of this parent's children's occupancy codes have been marked 
        %for pruning
        parents_with_all_children_toprune = parents((unique_count - child_count) == 0);
        if ~isempty(parents_with_all_children_toprune)
            %For all the parents_with_all_children_toprune, check if any of 
            %their children have previously been marked as leaf cells: if 
            %so, then we cannot prune the occupancy codes of the 
            %corresponding parents and we need to make ALL of the siblings 
            %of these children leaves as well
            parents_cannot_prune = [];
            for parent = parents_with_all_children_toprune'
                %Get all the children of the current parent
                children = ((myOT.FirstChildPtr{lvl - 1}(parent)):(myOT.FirstChildPtr{lvl - 1}(parent) + uint32(myOT.ChildCount{lvl - 1}(parent)) - 1));
                %Check if any of the children were marked as leaf cells
                if any(post_pruning_array{lvl}(children))
                    %Mark all the children of this parent as leaves
                    post_pruning_array{lvl}(children) = 1;
                    %We cannot prune the current parent's occupancy code, 
                    %so mark it to indicate that we will not be checking 
                    %this parent's ancestors 
                    parents_cannot_prune(end + 1) = parent;
                end
            end
            %Remove the parents whose occupancy codes cannot be pruned,
            %from the list
            parents_with_all_children_toprune(ismember(parents_with_all_children_toprune, parents_cannot_prune) == 1) = [];
        end %End check if ~isempty(parents_with_all_children_toprune)
        %For the remaining parents (if any), which have all of their 
        %children marked for pruning and none of these children are 
        %leaves ...
        if ~isempty(parents_with_all_children_toprune)
            %Mark these parents' occupancy codes for pruning
            toprune{lvl - 1}(((end + 1):(end + length(parents_with_all_children_toprune))), 1) = parents_with_all_children_toprune;
            %Also mark the occupancy codes of all the children of these 
            %parents for pruning from post_pruning_array, since none of 
            %these children will now be leaves (as they will be pruned 
            %away) 
            for parent = parents_with_all_children_toprune'
                %Get all the children of the current parent
                children = ((myOT.FirstChildPtr{lvl - 1}(parent)):(myOT.FirstChildPtr{lvl - 1}(parent) + uint32(myOT.ChildCount{lvl - 1}(parent)) - 1));
                %Mark these children's occupancy codes for pruning from
                %post_pruning_array
                toprune2{lvl}((end + 1):(end + length(children)), 1) = children;
                %Check if the current parent has all zero wavelet 
                %coefficients: if NOT, but all of its children do (since 
                %they were all marked for pruning and none of them are 
                %leaves), then the parent will be a leaf
                if ~any(all_zero_wav_cfs{lvl - 1} == parent)
                    post_pruning_array{lvl - 1}(parent) = 1; 
                end 
            end
        end %End check if ~isempty(parents_with_all_children_toprune)
        %Find the parents amongst "parents", whose children's occupancy 
        %codes have NOT all been marked for pruning
        parents_with_not_all_children_toprune = parents((unique_count - child_count) ~= 0);
        %For parents whose children have NOT all been marked for pruning 
        %(if any), we cannot prune the occupancy codes of these parents, so
        %all of the children of these parents, which have been marked for 
        %pruning (NOT the children that have not been marked for pruning) 
        %will become leaves
        if ~isempty(parents_with_not_all_children_toprune)
            for parent = parents_with_not_all_children_toprune'
                %Get all the children of the current parent
                children = ((myOT.FirstChildPtr{lvl - 1}(parent)):(myOT.FirstChildPtr{lvl - 1}(parent) + uint32(myOT.ChildCount{lvl - 1}(parent)) - 1));
                %Find the children that have been marked for pruning and 
                %set them to be leaves
                post_pruning_array{lvl}(children(ismember(children, toprune{lvl}) == 1)) = 1;
            end
        end %End check if ~isempty(parents_with_not_all_children_toprune)
    end %End check if ~isempty(toprune{lvl})
    disp(['Finished marking level ' num2str(lvl) ' for pruning']);
end %End lvl
disp('------------------------------------------------------------');

%Check if all the marked leaf cells have different control point signs on
%their corners; if not, make these cells internal and look for the next
%available leaf cells further down in the octree. We do not want leaf cells
%to have the same sign control points on all their corners, because then 
%these cells would not be subdivided at the decoder and their descendant 
%voxels would not be reconstructed.
disp('Checking if marked leaf cells have different control point signs on their corners ...');
%Find the first non-empty level in post_pruning_array
pp_first_nonempty = find(~cellfun(@isempty, post_pruning_array), 1);
%Find the first level in post_pruning_array, which contains leaf cells
for ppa_lvl = pp_first_nonempty:size(post_pruning_array, 1)
    test = find((post_pruning_array{ppa_lvl} == 1), 1);
    if ~isempty(test)
        ppa_first_leaf_level = ppa_lvl;
        break;
    end
end

for lvl = ppa_first_leaf_level:(max_lvl - 1)   %Assuming max_lvl = b + 1
    %Extract indices of all leaf cells at the current octree level
    leaves = find(post_pruning_array{lvl} == 1);
    %For each leaf cell, check its control point signs
    for leaf = leaves'
        current_ctrlpts = reconstructed_control_points{lvl}(ctrl_pts_pointers{lvl}((leaf*8 - 7):(leaf*8)));
        if abs(sum(sign(current_ctrlpts))) ~= 8
            %If the control point signs are not all the same, the cell can
            %remain a leaf
            continue;
        elseif (abs(sum(sign(current_ctrlpts))) == 8)||(~any(sign(current_ctrlpts)))
            %Otherwise, the cell cannot be a leaf, so set it to 0 in
            %post_pruning_array
            post_pruning_array{lvl}(leaf) = 0;
            %Remove this cell from toprune, because its children can no
            %longer be pruned away since they need to be checked for
            %potential leaves
            toprune{lvl}(toprune{lvl} == leaf) = [];
            %Check this cell's children (unless the children are at the 
            %voxel level) to see if they can be leaves (if their control 
            %point signs are different)
            if lvl < b
                current_children = (myOT.FirstChildPtr{lvl}(leaf)):(myOT.FirstChildPtr{lvl}(leaf) + uint32(myOT.ChildCount{lvl}(leaf)) - 1);
                children_not_for_leaves = [];
                for child = current_children
                    child_ctrlpts = reconstructed_control_points{lvl + 1}(ctrl_pts_pointers{lvl + 1}((child*8 - 7):(child*8)));
                    if (abs(sum(sign(child_ctrlpts))) == 8)||(~any(sign(child_ctrlpts)))
                        %If the control point signs are all the same, the 
                        %current child cannot be a leaf
                        children_not_for_leaves(end + 1) = child;
                        %Remove the child from toprune and toprune2, because
                        %its descendants will need to be checked for leaves and
                        %so its children can no longer be pruned away
                        toprune{lvl + 1}(toprune{lvl + 1} == child) = [];
                        toprune2{lvl + 1}(toprune2{lvl + 1} == child) = [];
                    elseif abs(sum(sign(child_ctrlpts))) ~= 8
                        %If the control point signs are different, the child
                        %can be a leaf, so set it to 1 in post_pruning_array
                        post_pruning_array{lvl + 1}(child) = 1;
                        %Remove the child from toprune2 as we do not want to
                        %prune away leaves
                        toprune2{lvl + 1}(toprune2{lvl + 1} == child) = [];   
                    end
                end %End child
                %If there are no children in children_not_for_leaves, no need
                %to check any more descendants of the current "leaf".
                %Otherwise, for all the children in children_not_for_leaves, 
                %check their descendants to find the next available leaves 
                %further down in the octree.
                if ~isempty(children_not_for_leaves)
                    check_other_descendants = children_not_for_leaves;
                    %For each level further down in the octree, check the
                    %descedants of the cells in check_other_descendants
                    for lvl2 = (lvl + 1):(b - 1)
                        descendants_not_leaves = [];
                        for desc = check_other_descendants
                            %Get the children of the current desc
                            desc_children = (myOT.FirstChildPtr{lvl2}(desc)):(myOT.FirstChildPtr{lvl2}(desc) + uint32(myOT.ChildCount{lvl2}(desc)) - 1);
                            %Check control points of each of the desc_children,
                            %to see if any of them can be leaves
                            for dc = desc_children
                                dc_ctrlpts = reconstructed_control_points{lvl2 + 1}(ctrl_pts_pointers{lvl2 + 1}((dc*8 - 7):(dc*8)));
                                if (abs(sum(sign(dc_ctrlpts))) == 8)||(~any(sign(dc_ctrlpts)))
                                    %If the control point signs are all the same, the 
                                    %current child cannot be a leaf
                                    descendants_not_leaves(end + 1) = dc;
                                    %Remove the child from toprune and toprune2, because
                                    %its descendants will need to be checked for leaves and
                                    %so its children can no longer be pruned away
                                    toprune{lvl2 + 1}(toprune{lvl2 + 1} == dc) = [];
                                    toprune2{lvl2 + 1}(toprune2{lvl2 + 1} == dc) = [];
                                elseif abs(sum(sign(dc_ctrlpts))) ~= 8
                                    %If the control point signs are different, the child
                                    %can be a leaf, so set it to 1 in post_pruning_array
                                    post_pruning_array{lvl2 + 1}(dc) = 1;
                                    %Remove the child from toprune2 as we do not want to
                                    %prune away leaves
                                    toprune2{lvl2 + 1}(toprune2{lvl2 + 1} == dc) = [];   
                                end
                            end %End dc
                        end %End desc
                        if ~isempty(descendants_not_leaves)
                            check_other_descendants = descendants_not_leaves;
                        else
                            break;
                        end
                    end %End lvl2 
                end %End check if ~isempty(children_not_for_leaves)
            end %End check if lvl < b
        end %End check if abs(sum(sign(current_ctrlpts))) ~= 8
    end %End leaf
end %End lvl
disp('------------------------------------------------------------');

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

%If there are any octree levels in post_pruning_array that contain no leaf
%cells, then clear the post_pruning_array at these levels because we do not
%need to transmit only 0 bits
for lvl = 1:size(post_pruning_array, 1)
    if ~isempty(post_pruning_array{lvl}) && (~any(post_pruning_array{lvl}))
        post_pruning_array{lvl} = [];
        disp(['Cleared level ' num2str(lvl) ' of post_pruning_array: no leaves']);
        disp('------------------------------------------------------------');
    end
end

ot_pruning_time = toc;
disp(' ');
disp('************************************************************');
disp(['Time taken to prune octree: ' num2str(ot_pruning_time) ' seconds']);
disp('************************************************************');


