%DOES WAVELET ANALYSIS FROM (LEAVES-1) TO THE ROOT LEVEL, WITHOUT ANY
%QUANTIZATION.

%Requires directory "phil": add this directory plus its sub-directories to
%the current MATLAB path.

ptcloud_file = '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10\boxer_voxelized10.ply';
b = 10;
vis_levels_ot = 4; %No. of octree levels for which we want to visualize the octree cell subdivision
vis_levels_ctrlpts = 2; %No. of octree levels for which we want to visualize the control point computation
vis_levels_wavelet = [];    %Octree level(s) for which we want to visualize the wavelet analysis steps

%Initialize a cell array to store the corner coordinates of each occupied
%cell at each level of the octree
corner_coords = cell(b, 1);
%Initialize a cell array to store only the unique corner coordinates at
%each octree level
unique_coords = cell(b, 1);
%Initialize a cell array to store the Bezier control points (signed
%distances) associated with the unique corner coordinates at each octree 
%level
control_points = cell(b, 1);
%Initialize a cell array of "pointers": there will be 8 pointers per
%occupied octree cell (1 per corner) at each octree level, which will point
%to the Bezier control point in control_points, associated with that
%corner. So, for each octree level "lvl", there will be 
%myOT.NodeCount(lvl) pointer arrays containining 8 elements each.  
ctrl_pts_pointers = cell(b, 1);

%-------------------------- Octree Construction --------------------------%

%Construct an octree of depth (b + 1), for the input point cloud
[myOT, mortonCodes_sorted, xyz_sorted] = construct_octree(ptcloud_file, b);

%Extract the set of occupied voxel coordinates (x, y, z) at all levels of
%the octree myOT 
[~, occupied_voxel_coords] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted);

%For each octree level (except the leaf level) ...
%for lvl = 1:b 
for lvl = 1:vis_levels_ot 
    disp('------------------------------------------------------------');
    disp(['Processing octree level ' num2str(lvl) ':']); 
    %Counter to keep track of how many corner coordinates we have stored at
    %the current octree level
    corner_coords_cntr = 1;
    disp('Computing corner coordinates for each occupied cell ...');
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Initialize a matrix to store the (x, y, z) coordinates of all the
        %corners of the current octree cell
        corners = zeros(8, 3);
        %Find the (x, y, z) coordinates of the origin of this cell (found
        %at the bottom left-hand corner farthest from the viewer)
        corners(1, :) = double(myOT.SpatialIndex{lvl}(occ_cell, :))*(2^(b + 1 - lvl)) - [0.5 0.5 0.5];
        %Find the (x, y, z) coordinates of the other 7 corners of this cell
        corners(2, :) = corners(1, :) + [2^(b + 1 - lvl) 0 0];
        corners(3, :) = corners(1, :) + [2^(b + 1 - lvl) 2^(b + 1 - lvl) 0];
        corners(4, :) = corners(1, :) + [0 2^(b + 1 - lvl) 0];
        corners(5, :) = corners(1, :) + [0 0 2^(b + 1 - lvl)];
        corners(6, :) = corners(1, :) + [2^(b + 1 - lvl) 0 2^(b + 1 - lvl)];
        corners(7, :) = corners(1, :) + [2^(b + 1 - lvl) 2^(b + 1 - lvl) 2^(b + 1 - lvl)];
        corners(8, :) = corners(1, :) + [0 2^(b + 1 - lvl) 2^(b + 1 - lvl)];
        %Store all the corner coordinates for the current cell in their
        %corresponding locations inside corner_coords
        corner_coords{lvl}(corner_coords_cntr:(corner_coords_cntr + 7), 1:3) = corners;   
        corner_coords_cntr = corner_coords_cntr + 8;
    end
    %Find all the unique corner coordinates at the current level of the
    %octree, and in ctrl_pts_pointers{lvl}, for each row that represents an
    %(x, y, z) coordinate from corner_coords{lvl}, store a pointer for this
    %coordinate into the unique_coords{lvl} array, to say what the
    %coordinate triplet in this row should be. These pointers will also act
    %as pointers to the control points stored in control_points{lvl}, since
    %control_points{lvl} will be the same size as unique_coords{lvl}, as
    %each unique corner coordinate triplet will have one corresponding 
    %Bezier control point (signed distance value) associated with it.
    [unique_coords{lvl}, ~, ctrl_pts_pointers{lvl}] = unique(corner_coords{lvl}, 'rows', 'stable');
end %end "lvl"

%--------------- Visualization of Octree Cell Subdivision ----------------%

disp('------------------------------------------------------------');
disp('Displaying octree levels ...'); 

%Read in and display the input point cloud (assume PLY format)
[~, ptcloud, ~] = plyRead(ptcloud_file);
figure;
scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 4)./255, ptcloud(:, 5)./255, ptcloud(:, 6)./255], 'filled');
axis equal; axis off;
hold on;

%Define a colour matrix for each octree level for which we wish to 
%visualize the octree cell subdivision
colours = hsv(vis_levels_ot);   %Each row contains 3 columns (R, G, B)

%Initialize a cell array to store, for each unique corner vertex at each
%octree level, a list of indices of the corner vertices that share an edge
%with that vertex. 
shared_edge_inds = cell(b, 1);

%Initialize an array to store handles for each octree level plot, to use
%for the legend 
h_ot = zeros(1, vis_levels_ot);

%Create a cell array of labels for the plotted octree levels, to use in the
%legend
otplot_legend = cell(1, vis_levels_ot);
for l = 1:vis_levels_ot
    if l == 1
        otplot_legend{l} = 'Level 1 (Root)';
    else
        otplot_legend{l} = ['Level ' num2str(l)];
    end
end

for lvl = 1:vis_levels_ot
    %Plot all the unique corner points at the current octree level, over
    %the top of the input 3D point cloud
    scatter3(unique_coords{lvl}(:, 1), unique_coords{lvl}(:, 2), unique_coords{lvl}(:, 3), 10, colours(lvl, :), 'filled');
    hold on;
    %For two of the unique vertices to share an edge, they must have at
    %least one coordinate (x, or y, or z) in common. Work out all the edges
    %for the unique vertices, and use them to display the horizonal and 
    %vertical octree grid lines that show the cell divisions at the current
    %octree level (i.e., connect the points that share each edge with a
    %straight line).
    for v = 1:size(unique_coords{lvl}, 1)
        %Initialize an array to store the indices of the shared vertices,
        %for the current vertex
        edge_pts_indices = [];
        %Initialize a counter to keep track of how many vertex indices
        %corresponding to shared edges have been found
        edge_pts_indices_cntr = 1;
        %Get the coordinates of the current vertex
        current_vtx = unique_coords{lvl}(v, :);
        %Find all the vertices (amongst the list of unique corner vertices
        %at the current octree level) that share an edge with current_vtx
        [x_same_rows, ~] = find(current_vtx(1) == unique_coords{lvl}(:, 1));
        x_same_rows(find(x_same_rows == v)) = [];   %Remove the current vertex index
        [y_same_rows, ~] = find(current_vtx(2) == unique_coords{lvl}(:, 2));
        y_same_rows(find(y_same_rows == v)) = [];   %Remove the current vertex index
        [z_same_rows, ~] = find(current_vtx(3) == unique_coords{lvl}(:, 3));
        z_same_rows(find(z_same_rows == v)) = [];   %Remove the current vertex index
        temp_cat = [x_same_rows; y_same_rows; z_same_rows];
        for v_ind = 1:length(temp_cat)
            %Only keep the vertex indices that appear twice in temp_cat
            if length(find(temp_cat == temp_cat(v_ind))) == 2
                edge_pts_indices(edge_pts_indices_cntr) = temp_cat(v_ind);
                edge_pts_indices_cntr = edge_pts_indices_cntr + 1;
            end
        end
        %Keep only the unique indices in edge_pts_indices
        edge_pts_indices = (unique(edge_pts_indices))';
        %Get the (x, y, z) coordinates corresponding to all the
        %edge_pts_indices
        edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
        %Connect the current_vtx to all the edge_pts_coords with straight
        %lines
        for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
            h_ot(lvl) = plot3([current_vtx(1) edge_pts_coords(nv, 1)], [current_vtx(2) edge_pts_coords(nv, 2)], [current_vtx(3) edge_pts_coords(nv, 3)], 'Color', colours(lvl, :));
            hold on;
        end
        %Store the current list of edge indices inside shared_edge_inds,
        %for future reference
        shared_edge_inds{lvl, v} = edge_pts_indices;
    end
end
hold off;
legend(h_ot, otplot_legend, 'Location', 'best');
title('Octree Subdivision and Computing Corner Coordinates of Octree Cells', 'Interpreter', 'none');

%-------------------- Computing Bezier Control Points --------------------%

%Initialize cell array to store the indices of octree cells that share each
%unique corner, at each octree level
shared_cells = cell(b, 1);
%Initialize a cell array to store the (x, y, z) coordinates of the nearest
%voxel found, for each unique corner coordinate, at each octree level
nearest_voxels = cell(b, 1);

for lvl = 1:b
    disp('------------------------------------------------------------');
    disp(['Computing Bezier control points for octree level ' num2str(lvl) ':']);
    %Counter to keep track of how many Bezier control points we have stored
    %at the current octree level
    control_points_cntr = 1;
    %For each unique corner coordinate, find the nearest occupied voxel in
    %any of the occupied octree cells at the current level that share this 
    %corner. Do this by measuring the Euclidean distance from the corner to 
    %the (centres of) voxels.
    disp('Computing signed distance to nearest occupied voxel, for each unique corner coordinate ...');
    for c = 1:size(unique_coords{lvl}, 1)
        %Find out which occupied octree cells share this corner
        row_inds = find(ismember(corner_coords{lvl}, unique_coords{lvl}(c, :), 'rows') == 1);
        shared_cells{lvl, c} = ceil(row_inds./8);        
        %Get the (x, y, z) coordinates of all the occupied voxels in the 
        %cells found above
        for i = 1:length(shared_cells{lvl, c})
            if i == 1
                current_occupied_voxel_coords = occupied_voxel_coords{lvl, shared_cells{lvl, c}(i)};
            else
                %Concatenate the coordinates from the cells that share the
                %current corner (concatenate along rows, so that we
                %maintain the 3-column structure)
                current_occupied_voxel_coords = cat(1, current_occupied_voxel_coords, occupied_voxel_coords{lvl, shared_cells{lvl, c}(i)});
            end
        end
        %Find the nearest neighbour to the current unique corner coordinate
        diff = unique_coords{lvl}(c, :) - current_occupied_voxel_coords; 
        euclid_dists = abs(arrayfun(@(idx) norm(diff(idx, :)), 1:size(diff, 1)));
        [min_euclid_dist, nearest_voxel_ind] = min(euclid_dists);
        %Store the (x, y, z) coordinates of the nearest voxel, for future
        %reference
        nearest_voxels{lvl}(c, 1:3) = current_occupied_voxel_coords(nearest_voxel_ind, :);
        %Check the direction of the normal for the nearest_voxel. If the 
        %normal is pointing away from the corresponding corner,
        %this indicates that the corner is "inside" the surface implied by the
        %point cloud in that octree cell, so make the distance value (computed
        %above) negative. If the normal is pointing towards the corner, this
        %implies that the corner is "outside" the surface, so leave the
        %distance value positive. Store the signed distance value in 
        %control_points: this represents the Bezier control point associated 
        %with the corresponding cell corner.
%         if <normal is pointing away)
%             control_points{lvl}(control_points_cntr) = -min_euclid_dist;
%         elseif <normal is pointing towards>
%             control_points{lvl}(control_points_cntr) = min_euclid_dist;
%         end
        control_points{lvl}(control_points_cntr) = min_euclid_dist; %****TEMPORARY CODE ONLY, UNTIL I GET NORMALS****
        control_points_cntr = control_points_cntr + 1;   
    end %End unique corner coordinate 
    %Arrange all the control points at the current octree level, into a 
    %column vector instead of a row vector (purely for visualization 
    %reasons: it is easier to scroll through a long column vector than a 
    %row vector)
    control_points{lvl} = (control_points{lvl})';
end %End octree level
    
%------------- Visualization of Control Point Computation ----------------%

% disp('------------------------------------------------------------');
% disp('Visualizing control point computation ...'); 
% 
% %for lvl = 1:vis_levels_ctrlpts   
% for lvl = 3   
%     %For each unique corner coordinate ...
%     for c = 1:size(unique_coords{lvl}, 1)
%         %Display the input point cloud
%         figure;
%         scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 4)./255, ptcloud(:, 5)./255, ptcloud(:, 6)./255], 'filled');
%         axis equal; axis off;
%         hold on; 
%         %Display all the octree cells at this level
%         for j = 1:size(unique_coords{lvl}, 1)
%             %Get the list of vertex indices that share an edge with the current
%             %vertex
%             edge_pts_indices = shared_edge_inds{lvl, j};
%             %Get the (x, y, z) coordinates corresponding to all the
%             %edge_pts_indices
%             edge_pts_coords = unique_coords{lvl}(edge_pts_indices, :);
%             %Connect the current vertex to all the edge_pts_coords with 
%             %straight lines
%             for nv = 1:length(edge_pts_indices) %"nv" stands for neighbouring vertex
%                 plot3([unique_coords{lvl}(j, 1) edge_pts_coords(nv, 1)], [unique_coords{lvl}(j, 2) edge_pts_coords(nv, 2)], [unique_coords{lvl}(j, 3) edge_pts_coords(nv, 3)], 'Color', colours(lvl, :));
%                 hold on;
%             end
%         end
%         %Initialize a matrix to store all the corner coordinates of all the 
%         %octree cells at the current level, which share the current vertex 
%         shared_cell_coords = [];
%         %Get the indices of the shared octree cells (at the current level)
%         %for the current vertex
%         shared_cell_inds = shared_cells{lvl, c};  
%         %Get all the corner coordinates for all the shared cells
%         for sc = 1:length(shared_cell_inds)
%             if sc == 1
%                 shared_cell_coords = corner_coords{lvl}((shared_cell_inds(sc)*8 - 7):(shared_cell_inds(sc)*8), :);
%             else
%                 shared_cell_coords = cat(1, shared_cell_coords, corner_coords{lvl}((shared_cell_inds(sc)*8 - 7):(shared_cell_inds(sc)*8), :));
%             end
%         end
%         %Initialize a cell array to store the indices of the edge vertices
%         %for the shared cells of the current vertex
%         sc_edge_pts_indices = cell(size(shared_cell_coords, 1), 1);
%         %For each corner of the shared cells, figure out which other
%         %vertices in shared_cell_coords it is connected to
%         for scc = 1:size(shared_cell_coords, 1)
%             %Initialize a counter to keep track of how many shared corner
%             %coordinates have been found so far
%             sc_edge_pts_cntr = 1;
%             [x_same_rows, ~] = find(shared_cell_coords(scc, 1) == shared_cell_coords(:, 1));
%             [y_same_rows, ~] = find(shared_cell_coords(scc, 2) == shared_cell_coords(:, 2));
%             [z_same_rows, ~] = find(shared_cell_coords(scc, 3) == shared_cell_coords(:, 3));  
%             temp_cat = [x_same_rows; y_same_rows; z_same_rows];
%             %Only keep the vertex indices that appear twice in temp_cat
%             for v_ind = 1:length(temp_cat)
%                 if length(find(temp_cat == temp_cat(v_ind))) == 2
%                     sc_edge_pts_indices{scc}(sc_edge_pts_cntr) = temp_cat(v_ind);
%                     sc_edge_pts_cntr = sc_edge_pts_cntr + 1;
%                 end
%             end
%             %Keep only the unique indices in sc_edge_pts_indices{scc}
%             sc_edge_pts_indices{scc} = (unique(sc_edge_pts_indices{scc}))';
%         end
%         %Connect all the corner vertices of the shared octree cells, to
%         %other corner vertices on these cells, with which they share an 
%         %edge. The result will be a thick, black outline of the octree 
%         %cells at the current octree level, which share the current corner 
%         %vertex.
%         for c1 = 1:size(sc_edge_pts_indices, 1)
%             for c2 = 1:numel(sc_edge_pts_indices{c1})
%                 %Get the (x, y, z) coordinates of current two vertices that
%                 %will be connected
%                 xyz_c1 = shared_cell_coords(c1, :);
%                 xyz_c2 = shared_cell_coords(sc_edge_pts_indices{c1}(c2), :);
%                 h(2) = plot3([xyz_c1(1) xyz_c2(1)], [xyz_c1(2) xyz_c2(2)], [xyz_c1(3) xyz_c2(3)], 'k', 'LineWidth', 3);
%                 hold on;
%             end
%         end
%         title({['OCTREE LEVEL ' num2str(lvl) ':'], 'Identifying Shared Octree Cells and', 'Computing Distance From Each Corner to the Nearest Voxel in Those Cells'}, 'Interpreter', 'none');
%         %Plot the current corner point (in red, with a black outline) on 
%         %top of the input 3D point cloud and the octree grid
%         h(1) = scatter3(unique_coords{lvl}(c, 1), unique_coords{lvl}(c, 2), unique_coords{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
%         hold on;
%         %Colour in the nearest voxel (in blue, with a black outline) found
%         %for the current corner point. Note that the size of this voxel 
%         %will be exaggerated on this plot, in order to make it easier to 
%         %see.  
%         h(3) = scatter3(nearest_voxels{lvl}(c, 1), nearest_voxels{lvl}(c, 2), nearest_voxels{lvl}(c, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
%         legend(h, 'Current corner point', 'Shared octree cells (at current octree level)', 'Nearest voxel point', 'Location', 'best');
%         hold off;
%     end %End current unique corner at level "lvl"
% end %End octree level "lvl"
% disp('');


%------------------------- Wavelet Decomposition -------------------------%

%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) across all 
%octree blocks and levels
wavelet_coeffs = cell(b, 1);
%Initialize a flag that indicates whether or not we want to display
%visualizations at a given octree level (0 => no; 1 => yes)
w_vis_flag = 0;

%For each octree level, starting from one level before the leaves and
%working up to the root ...
%for lvl = b:-1:2
for lvl = 4:-1:2
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    cnr_coords_cntr = 1;
    %Initialize a counter to keep track of the number of wavelet
    %coefficients produced at this level
    wavelet_coeffs_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Extract the cell's 8 corner coordinates
        cell_corner_coords = corner_coords{lvl}(cnr_coords_cntr:(cnr_coords_cntr + 7), :);
        %Get the pointer to the parent cell of the current octree cell
        %(this parent will be found at the previous octree level)
        parent_ptr = myOT.ParentPtr{lvl}(occ_cell);
        %Get the 8 corner coordinates of the parent cell 
        parent_corner_coords = corner_coords{lvl - 1}((parent_ptr*8 - 7):parent_ptr*8, :); 
        %For each corner of the current cell ...
        for cnr = 1:8
            %Check if the current corner already has a wavelet coefficient 
            %associated with it (since some corners will be shared amongst 
            %different octree cells at the same level) 
            possible_inds = ctrl_pts_pointers{lvl}(cnr_coords_cntr:(cnr_coords_cntr + 7));
            if ~isempty(wavelet_coeffs{lvl})
                if length(wavelet_coeffs{lvl}) >= possible_inds(cnr) 
                    if ~isempty(wavelet_coeffs{lvl}(possible_inds(cnr)))
                        continue;          
                    end
                end
            end
            %Flag to indicate if current corner is on a parent edge (0 =>
            %no; 1 => yes)
            on_p_edge = 0;
            %Flag to indicate if current corner is on a parent face (0 =>
            %no; 1 => yes)
            on_p_face = 0;
            %Display visualizations only for the chosen octree level(s)
            if ~isempty(find(lvl == vis_levels_wavelet))
                w_vis_flag = 1;
                figure;
                %Plot the input point cloud
                scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 4)./255, ptcloud(:, 5)./255, ptcloud(:, 6)./255], 'filled');
                axis equal; axis off;
                hold on;
                %Outline the current octree cell in thick black lines 
                %(i.e., connect all the vertices in cell_corner_coords to 
                %the other vertices in cell_corner_coords,with which they 
                %share an edge). Since we are only plotting one cell, we 
                %know in advance how the vertices are connected, so first 
                %construct a matrix of edge indices.
                edges_onecell = zeros(8, 3);    %Each corner vertex is connected to 3 others
                edges_onecell(1, :) = [2, 4, 5];    %Corner 1 is connected to corners 2, 4, and 5
                edges_onecell(2, :) = [1, 3, 6];
                edges_onecell(3, :) = [2, 4, 7];
                edges_onecell(4, :) = [1, 3, 8];
                edges_onecell(5, :) = [1, 6, 8];
                edges_onecell(6, :) = [2, 5, 7];
                edges_onecell(7, :) = [3, 6, 8];
                edges_onecell(8, :) = [4, 5, 7];
                %Plot the root octree cell, for reference
                for rootv1 = 1:8
                    for rootv2 = 1:3
                        h_w(1) = plot3([corner_coords{1}(rootv1, 1) corner_coords{1}(edges_onecell(rootv1, rootv2), 1)], [corner_coords{1}(rootv1, 2) corner_coords{1}(edges_onecell(rootv1, rootv2), 2)], [corner_coords{1}(rootv1, 3) corner_coords{1}(edges_onecell(rootv1, rootv2), 3)], 'r');
                        hold on;
                    end
                end
                %Connect the vertices in the current cell, according to
                %edges_onecell
                for cvtx1 = 1:8
                    for cvtx2 = 1:3
                        h_w(2) = plot3([cell_corner_coords(cvtx1, 1) cell_corner_coords(edges_onecell(cvtx1, cvtx2), 1)], [cell_corner_coords(cvtx1, 2) cell_corner_coords(edges_onecell(cvtx1, cvtx2), 2)], [cell_corner_coords(cvtx1, 3) cell_corner_coords(edges_onecell(cvtx1, cvtx2), 3)], 'k', 'LineWidth', 3);
                        hold on;
                    end
                end
                %Outline the parent cell in thick, dotted black lines, by 
                %connecting the vertices in the parent cell according to
                %edges_onecell
                for pvtx1 = 1:8
                    for pvtx2 = 1:3
                        h_w(3) = plot3([parent_corner_coords(pvtx1, 1) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 1)], [parent_corner_coords(pvtx1, 2) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 2)], [parent_corner_coords(pvtx1, 3) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 3)], '--k', 'LineWidth', 3);
                        hold on;
                    end
                end               
                %Plot the current corner point (in red, with a black 
                %outline) on the current figure
                h_w(4) = scatter3(cell_corner_coords(cnr, 1), cell_corner_coords(cnr, 2), cell_corner_coords(cnr, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
                hold on;   
            end %End first part of the visualization code for the current corner
            
            %For the current corner, check if it is on a parent edge (i.e.,
            %if it shares at least 2 same coordinates (out of x, y, or z)
            %with two of the parent vertices), or on a parent face (i.e.,
            %only 1 of its coordinates (usually y, but depends on input
            %coordinate system) is the same as four parents' coordinates 
            %(must have same coordinate in common: either x, or y, or z)  
            parent_row_inds = [];
            [x_same, ~] = find(cell_corner_coords(cnr, 1) == parent_corner_coords(:, 1));
            [y_same, ~] = find(cell_corner_coords(cnr, 2) == parent_corner_coords(:, 2));
            [z_same, ~] = find(cell_corner_coords(cnr, 3) == parent_corner_coords(:, 3));
            temp_cat = [x_same; y_same; z_same];
            if ((isempty(x_same) + isempty(y_same) + isempty(z_same)) == 2)
                %Corner is on a parent face
                on_p_face = 1;
                %Get the row indices of the parent vertices on this face
                if ~isempty(x_same)
                    parent_row_inds = x_same;
                elseif ~isempty(y_same)
                    parent_row_inds = y_same;
                elseif ~isempty(z_same)
                    parent_row_inds = z_same;
                end
            else
                %Check if any vertex indices appear twice in temp_cat
                temp_cntr = 1;
                for v = 1:length(temp_cat)
                    if length(find(temp_cat == temp_cat(v))) == 2
                        parent_row_inds(temp_cntr) = temp_cat(v);
                        temp_cntr = temp_cntr + 1;
                    end
                end
                if (length(unique(parent_row_inds)) == 2)
                    %Corner is on a parent edge
                    on_p_edge = 1;
                    %Keep only the unique indices in parent_row_inds
                    parent_row_inds = unique(parent_row_inds);
                end      
            end
            
            %If this corner's coordinates are the same as one of the
            %corner coordinates of the parent
            if sum(ismember(parent_corner_coords, cell_corner_coords(cnr, :), 'rows') > 0)
                %Do nothing (the signal on this vertex is a low-pass
                %coefficient), but still leave an empty place (MATLAB will
                %automatically insert a 0 here) in the wavelet_coeffs cell 
                %array, so that the wavelet coefficients in there 
                %correspond to the same locations in unique_coords
                wavelet_coeffs_cntr = wavelet_coeffs_cntr + 1;
                continue;
            %If this corner vertex lies on a parent edge
            elseif on_p_edge == 1
                %Get the indices of the rows of these coordinates in
                %parent_corner_coords (should be only 2 rows) 
                %parent_row_inds = find(sum(ismember(parent_corner_coords, cell_corner_coords(cnr, :)), 2) == 2);
                if w_vis_flag == 1
                    for pc = 1:length(parent_row_inds)
                        %Circle the parent corners (in blue)
                        h_w(5) = scatter3(parent_corner_coords(parent_row_inds(pc), 1), parent_corner_coords(parent_row_inds(pc), 2), parent_corner_coords(parent_row_inds(pc), 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                        hold on;
                    end
                end
                %Get the Bezier control points stored at the corner, 
                %vertices defined by parent_row_inds
                all_ctrlpts_ptrs = ctrl_pts_pointers{lvl-1}((parent_ptr*8 - 7):parent_ptr*8, :);
                ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                ctrlpt1 = control_points{lvl-1}(ctrlpt1_ptr);
                ctrlpt2 = control_points{lvl-1}(ctrlpt2_ptr);
                %Average the signal (Bezier control points) on the 2
                %vertices of the parent edge
                avg_signal = (ctrlpt1 + ctrlpt2)/2;
            %If this corner vertex lies on a parent face 
            elseif on_p_face == 1
                %Get the indices of the rows of these coordinates in
                %parent_corner_coords (should be 4 rows) 
                %parent_row_inds = find(sum(ismember(parent_corner_coords, cell_corner_coords(cnr, :)), 2) == 1);
                if w_vis_flag == 1
                    for pc = 1:length(parent_row_inds)
                        %Circle the parent corners (in blue)
                        h_w(5) = scatter3(parent_corner_coords(parent_row_inds(pc), 1), parent_corner_coords(parent_row_inds(pc), 2), parent_corner_coords(parent_row_inds(pc), 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                        hold on;
                    end
                end
                %Get the Bezier control points stored at the corner
                %vertices defined by parent_row_inds
                all_ctrlpts_ptrs = ctrl_pts_pointers{lvl-1}((parent_ptr*8 - 7):parent_ptr*8, :);
                ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                ctrlpt3_ptr = all_ctrlpts_ptrs(parent_row_inds(3));
                ctrlpt4_ptr = all_ctrlpts_ptrs(parent_row_inds(4));
                ctrlpt1 = control_points{lvl-1}(ctrlpt1_ptr);
                ctrlpt2 = control_points{lvl-1}(ctrlpt2_ptr);
                ctrlpt3 = control_points{lvl-1}(ctrlpt3_ptr);
                ctrlpt4 = control_points{lvl-1}(ctrlpt4_ptr);
                %Average the signal (Bezier control points) on the 4 
                %vertices of the parent face
                avg_signal = (ctrlpt1 + ctrlpt2 + ctrlpt3 + ctrlpt4)/4;  
            %If this corner vertex lies somewhere in the centre of the
            %parent's block (i.e., neither on a parent's edge or on a
            %parent's face)
            else                
                if w_vis_flag == 1
                    for pc = 1:8
                        %Circle the parent corners (in blue)
                        h_w(5) = scatter3(parent_corner_coords(pc, 1), parent_corner_coords(pc, 2), parent_corner_coords(pc, 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                        hold on;
                    end
                end
                %Get the Bezier control points stored at all of the corner
                %vertices (8) of the parent cell
                ctrlpts_pointers = ctrl_pts_pointers{lvl-1}((parent_ptr*8 - 7):parent_ptr*8, :);
                ctrlpts = control_points{lvl-1}(ctrlpts_pointers);
                %Average the signal (Bezier control points) on the 8
                %vertices of the parent block
                avg_signal = mean(ctrlpts);
            end
            %Get the Bezier control point stored at the current corner of
            %the current octree cell (the child cell)
            child_ctrlpt = ctrl_pts_pointers{lvl}(occ_cell*8 - 8 + cnr);
            %Subtract the average from the signal (Bezier control point) at
            %the current child vertex. The result is the transform 
            %(wavelet) coefficient of the child vertex.
            wavelet_coeffs{lvl}(wavelet_coeffs_cntr) = child_ctrlpt - avg_signal;
            wavelet_coeffs_cntr = wavelet_coeffs_cntr + 1;
            if w_vis_flag == 1
                hold off;
                legend(h_w, 'Root cell (for reference)', 'Current cell', 'Parent cell', 'Current corner', 'Parent corners', 'Location', 'best');
                title({'Computing Wavelet Coefficients', ['for Each Occupied Octree Cell at Level ' num2str(lvl)]});
                w_vis_flag = 0; %Reset flag for next octree cell
            end
        end %End current corner (of current cell)
        cnr_coords_cntr = cnr_coords_cntr + 8;
    end %End current octree cell
    %Arrange all the wavelet coefficients at the current octree level, into 
    %a column vector instead of a row vector (purely for visualization 
    %reasons: it is easier to scroll through a long column vector than a 
    %row vector). Also, wavelet_coeffs will then be in the same format as
    %the control_points cell array.
    wavelet_coeffs{lvl} = (wavelet_coeffs{lvl})';
end %End current octree level


%---------------------------- Reconstruction -----------------------------%









