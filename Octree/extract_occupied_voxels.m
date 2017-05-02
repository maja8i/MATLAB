%-------------------------------------------------------------------------%

%Extract the occupied voxel coordinates, and their corresponding normals
%(if the input point cloud has normals) in each occupied cell at each level
%of an octree, and plot the occupied voxels from each cell in a different
%colour on the same plot.

%---- INPUTS ----

%myOT: Octree class, with a number of different properties (can be obtained
%      by using construct_octree()).

%mortonCodes_sorted: Morton codes representing the (x, y, z) locations in
%                    the input point cloud, sorted in ascending order (can
%                    be obtained from construct_octree()).

%xyz_sorted: Input (x, y, z) triplets arranged in the same order as their
%            corresponding Morton codes (can be obtained from
%            construct_octree()).

%Optional inputs:

%normals_sorted: Input normals (also (x, y, z) triplets) arranged in the
%                same order as their corresponding (x, y, z) location
%                values.

%---- OUTPUTS ----

%occupied_Morton_codes: Cell array containing Morton codes for the occupied
%                       voxels in each occupied cell, at each level of the 
%                       octree.

%occupied_voxel_coords: Cell array containing the x, y, z coordinates of 
%                       the occupied voxels in each occupied cell, at each
%                       level of the octree.

%occupied_voxel_normals: Cell array containing the x, y, z coordinates of
%                        the normals of each of the occupied voxels in each
%                        occupied cell, at each level of the octree.

%-------------------------------------------------------------------------%

function [occupied_Morton_codes, occupied_voxel_coords, occupied_voxel_normals] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted, varargin)

%Initialize cell arrays to store Morton codes and corresponding x, y, z
%coordinates for occupied voxels in each occupied cell at each level of the
%octree
occupied_Morton_codes = cell((myOT.Depth + 1), max(myOT.NodeCount));
occupied_voxel_coords = cell((myOT.Depth + 1), max(myOT.NodeCount));

%If the input point cloud has normals and these have been passed to the
%current function
if numel(varargin) == 1
    normals_sorted = varargin{1};
    %Initialize a cell array to store the normals for all of the occupied
    %voxels in each occupied cell at each level of the octree, in the same
    %order as their corresponding voxel coordinates and Morton codes in the 
    %occupied_voxel_coords and occupied_Morton_codes cell arrays above
    occupied_voxel_normals = cell((myOT.Depth + 1), max(myOT.NodeCount));
end

%For each level of the octree myOT ...
for lvl = 1:(myOT.Depth + 1)  
    if numel(varargin) == 1
        disp(['Extracting occupied voxel coordinates and associated normal vectors in all occupied cells at octree level ' num2str(lvl) ' ...']);
    else
         disp(['Extracting occupied voxel coordinates in all occupied cells at octree level ' num2str(lvl) ' ...']);
    end
    disp(['Total no. of occupied cells to process: ' num2str(myOT.NodeCount(lvl))]);
    disp('------------------------------------------------------------');
    %For each occupied cell at this level ...
    for n = 1:myOT.NodeCount(lvl)   %NodeCount gives the no. of occupied cells at each octree level   
        %Get the index of the first occupied voxel in the n-th occupied 
        %cell at the current level ("lvl")
        begin = myOT.FirstDescendantPtr{lvl}(n);   
        %Get the total number of occupied voxels in the n-th occupied cell
        %at the current level ("lvl")
        count = myOT.DescendantCount{lvl}(n);
        %Get the list of Morton codes of all the occupied voxels in the
        %n-th occupied cell at the current level ("lvl")
        current_Morton_set = mortonCodes_sorted(begin:(begin + count - 1));  
        occupied_Morton_codes{lvl, n} = current_Morton_set;
        %Get the corresponding x, y, z location values for the Morton codes
        %in current_Morton_set
        current_voxel_set = xyz_sorted(begin:(begin + count - 1), :);
        occupied_voxel_coords{lvl, n} = current_voxel_set;
        if numel(varargin) == 1
            %Get the corresponding normal x, y, z values for the occupied
            %voxels found above
            occupied_voxel_normals{lvl, n} = normals_sorted(begin:(begin + count - 1), :);
        end
    end
end

% %For each level of the octree myOT ...
% for lvl = 1:(myOT.Depth + 1)
%     %Define a matrix of colours, one colour per occupied cell at the 
%     %current octree level
%     cell_colours = hsv(myOT.NodeCount(lvl));   %Each row contains 3 columns (R, G, B)
%     %Plot the occupied voxels in each occupied cell, in a different colour
%     %as defined in cell_colours
%     figure;
%     disp('------------------------------------------------------------');
%     disp(['Plotting occupied cell voxels at octree level ' num2str(lvl) ' ...']);
%     for n = 1:myOT.NodeCount(lvl)
%         plot3(occupied_voxel_coords{lvl, n}(:, 1), occupied_voxel_coords{lvl, n}(:, 2), occupied_voxel_coords{lvl, n}(:, 3), '.', 'Color', cell_colours(n, :), 'MarkerSize', 5);
%         axis equal;
%         if lvl == 1
%             title(['Octree Level ' num2str(lvl) ' (Root)']);
%         else
%             title(['Octree Level ' num2str(lvl)]);
%         end
%         hold on;
%     end
% end