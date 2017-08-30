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

%The below code prunes branches (occupancy codes) of octree cells that 
%have all zero wavelet coefficients throughout the branch, so the 
%pruning is done bottom-up from the voxel level ...
tic;
for lvl = max_lvl %Here assume max_lvl = b + 1
    %If there exist any occupied voxels, which contain all zero wavelet 
    %coefficients ...
    if ~isempty(all_zero_wav_cfs{lvl})
        %disp('Found voxels with all-zero wavelet coefficients: ');
        %disp(num2str(all_zero_wav_cfs{lvl}'));
        %For each of these voxels ...
        for az_cell = all_zero_wav_cfs{lvl}
            %Find the parent of the current voxel
            parent_cell = myOT.ParentPtr{lvl}(az_cell);
            %Find all the children (voxels) of parent_cell (one of these
            %children will be the current voxel)
            children = ((myOT.FirstChildPtr{lvl - 1}(parent_cell)):(myOT.FirstChildPtr{lvl - 1}(parent_cell) + uint32((myOT.ChildCount{lvl - 1}(parent_cell))) - 1));
            %Check if any of the other children (voxels) also have all zero 
            %wavelet coefficients
            az_children = children(find(ismember(children, all_zero_wav_cfs{lvl}) == 1));
            %If ALL of the children have all zero wavelet coefficients
            if length(az_children) == length(children)
                %Mark the occupancy code of parent_cell for pruning
                toprune{lvl - 1}(length(toprune{lvl - 1}) + 1) = parent_cell;
                %If the parent_cell also has all zero wavelet coefficients
                if (ismember(parent_cell, all_zero_wav_cfs{lvl - 1}) == 1)
                    %Successively check ancestors of the current
                    %parent_cell, to see if we can prune further up the 
                    %octree 
                    for lvl2 = (lvl - 1):-1:(start_lvl + 1)
                        %Get the parent of the current parent_cell
                        parent2 = myOT.ParentPtr{lvl2}(parent_cell);
                        %Get the children of the grandparent (parent2) -
                        %one of these children will be the current
                        %parent_cell
                        children2 = ((myOT.FirstChildPtr{lvl2 - 1}(parent2)):(myOT.FirstChildPtr{lvl2 - 1}(parent2) + uint32((myOT.ChildCount{lvl2 - 1}(parent2))) - 1));
                        %Check if any az_children2 have all zero wavelet
                        %coefficients
                        az_children2 = children2(find(ismember(children2, all_zero_wav_cfs{lvl2}) == 1));
                        %If ALL of the children have all zero wavelet
                        %coefficients
                        if(length(az_children2) == length(children2))
                            %Mark the occupancy code for the current 
                            %ancestor parent for pruning
                            toprune{lvl2 - 1}(length(toprune{lvl2 - 1}) + 1) = parent2;
                            %Mark the current parent_cell for pruning (from
                            %post_pruning_array), since it will not be a
                            %leaf
                            toprune2{lvl2}(length(toprune2{lvl2}) + 1) = parent_cell;
                            %Mark the occupancy codes of parent_cell's
                            %siblings (i.e., parent2's other occupied 
                            %children) for pruning too, if these cells 
                            %have occupancy codes (occupied children) of
                            %their own
                            if ~isempty(myOT.OccupancyCode{lvl2}(children2))
                                toprune{lvl2}((length(toprune{lvl2}) + 1):(length(toprune{lvl2}) + length(children2))) = children2;
                                %These children will not be leaves either
                                %(since they will be pruned off), so mark 
                                %them for pruning from post_pruning_array
                                toprune2{lvl2}((length(toprune2{lvl2}) + 1):(length(toprune2{lvl2}) + length(children2))) = children2;
                            end
                            %If parent2 also has all zero wavelet 
                            %coefficients (as do all of its children) 
                            if (ismember(parent2, all_zero_wav_cfs{lvl2 - 1}) == 1)
                                %Continue checking other ancestors, with 
                                %the current parent2 being the new 
                                %parent_cell
                                parent_cell = parent2;
                                continue;  
                            %If the ancestor parent does NOT have all zero
                            %wavelet coefficients (but all of its children
                            %do)
                            else
                                %The ancestor parent will become a leaf
                                post_pruning_array{lvl2 - 1}(parent2) = 1;
                                %Mark the current parent_cell for
                                %pruning from post_pruning_array, since
                                %it is not a leaf
                                toprune2{lvl2}(length(toprune2{lvl2}) + 1) = parent_cell;    
                            end
                        %If NOT all of parent2's children (parent_cell's 
                        %siblings) have all zero wavelet coefficients
                        else
                            %We cannot prune the occupancy code of parent2,
                            %but the current parent_cell will become a leaf
                            post_pruning_array{lvl2}(parent_cell) = 1;
                            %Stop searching other ancestors of the current
                            %voxel; move on to the next voxel (az_cell)
                            break;
                        end 
                    end %End lvl2
                %If parent_cell does NOT have all zero wavelet coefficients 
                %(but all of its children do)    
                else
                    %The current parent_cell will become a leaf after
                    %pruning its children away
                    post_pruning_array{lvl - 1}(parent_cell) = 1;
                    %Do not need to check ancestors of parent_cell, so move
                    %on to checking the next voxel (az_cell) instead
                    continue;
                end
            %If NOT ALL of the children (voxels) of the current parent_cell
            %have all zero wavelet coefficients
            else
                %We cannot prune the occupancy code of the current 
                %parent_cell, so move on to checking the next voxel 
                %(az_cell) and its parent
                continue;
            end          
        end %End az_cell
    end %End check if ~isempty(all_zero_wav_cfs{lvl})
end %End lvl    

%For each octree level covered in toprune ...
for lvl = 1:size(toprune, 1)
    if ~isempty(toprune{lvl})
        %Keep only the unique indices
        toprune{lvl} = unique(toprune{lvl}, 'stable');
        %Make toprune{lvl} a column vector rather than a row vector (only 
        %because column vectors are easier to scroll through)
        toprune{lvl} = toprune{lvl}';
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
        %Keep only the unique indices 
        toprune2{lvl} = unique(toprune2{lvl}, 'stable');
        %Make toprune2{lvl} a column vector rather than a row vector (only
        %because column vectors are easier to scroll through)
        toprune2{lvl} = toprune2{lvl}';
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