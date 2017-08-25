function [wavelet_coeffs, reconstructed_control_points] = wavelet_analysis_loop(myOT, corner_coords, control_points, ctrl_pts_pointers, start_lvl, max_lvl, b, q_stepsize, ptcloud_file)

%Initialize a cell array to store the transform (wavelet) coefficients for
%all the unique corner vertices (1 coefficient per vertex) across all 
%octree blocks and levels, starting from start_lvl and going up to one
%level before the leaves
wavelet_coeffs = cell((b + 1), 1);  %These will be quantized coefficients
%Initialize a cell array to store the reconstructed signal at each vertex
%of each octree cell at every level from start_lvl to one level before the
%leaves
reconstructed_control_points = cell(size(control_points, 1), 1);
%Initialize a flag that indicates whether or not we want to display
%visualizations at a given octree level (0 => no; 1 => yes). IMPORTANT:
%This must always be initialized to 0; it will be automatically set to 1
%later, if the vis_levels_wavelet array (see below) is not empty.
w_vis_flag = 0;
%Octree level(s) for which we want to visualize the wavelet analysis steps
%(leave the below array empty if you do not wish to visualize any wavelet 
%analysis steps) 
vis_levels_wavelet = [];    

if ~isempty(vis_levels_wavelet)
    %Read in the input point cloud (assume PLY format)
    [~, ptcloud, ~] = plyRead(ptcloud_file);
end

%Quantize all the control points at octree level start_lvl
for cpt = 1:size(control_points{start_lvl}, 1)
    control_points{start_lvl}(cpt, :) = quantize_uniform_scalar(control_points{start_lvl}(cpt, :), q_stepsize);
end
%Add these control points to reconstructed_control_points
reconstructed_control_points{start_lvl} = control_points{start_lvl};

%For each octree level, starting from start_lvl and working up to 1 level
%before the leaves ...
%tic;
%for lvl = start_lvl:b
for lvl = start_lvl:(max_lvl - 1)
    disp(['Computing wavelet coefficients between octree levels ' num2str(lvl) ' and ' num2str(lvl + 1) ' ...']);
    disp('------------------------------------------------------------');
    tic;
    %Initialize a counter for the corner coordinates of the occupied cells 
    %at this level 
    parent_cnr_coords_cntr = 1;
    %For each occupied octree cell at the current level ...
    for occ_cell = 1:myOT.NodeCount(lvl)
        %Extract the cell's 8 corner coordinates. This cell will represent
        %our parent cell at the current level, since we will compute
        %wavelet coefficients for its child cells.
        parent_corner_coords = corner_coords{lvl}(parent_cnr_coords_cntr:(parent_cnr_coords_cntr + 7), :);
        %Get the pointer to the first child cell of the current parent cell
        child_ptr = myOT.FirstChildPtr{lvl}(occ_cell);
        %For each child cell of the current cell ...
        for child = 1:myOT.ChildCount{lvl}(occ_cell)
            %Get the 8 corner coordinates of the current child cell. NOTE:
            %Below, "child" is converted to type double, because it is
            %uint8 by default, and if the result of any of the 
            %multiplications is > 255, the answer will be truncated and the
            %final result will be incorrect.
            child_corner_coords = corner_coords{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
            %For each corner of the current child cell ...
            for cnr = 1:8
                %Check if the current corner already has a wavelet 
                %coefficient associated with it (since some corners will be 
                %shared amongst different octree cells and we want to make
                %sure that we process only the UNIQUE corner vertices at 
                %each octree level) 
                possible_inds = ctrl_pts_pointers{lvl + 1}(((child_ptr - 1)*8 + 8*double(child) - 7):(child_ptr*8 + 8*double(child) - 8), :);
                if ~isempty(wavelet_coeffs{lvl + 1})
                    if length(wavelet_coeffs{lvl + 1}) >= possible_inds(cnr) 
                        if ~isempty(wavelet_coeffs{lvl + 1}(possible_inds(cnr)))
                        %if (wavelet_coeffs{lvl + 1}(possible_inds(cnr))) ~= 0
                            continue;          
                        end
                    end
                end  
                %Flag to indicate if current corner is on a parent edge 
                %(0 => no; 1 => yes)
                on_p_edge = 0;
                %Flag to indicate if current corner is on a parent face 
                %(0 => no; 1 => yes)
                on_p_face = 0;
                
                %--------------------- Visualization ---------------------%
                
                %Display visualizations only for the chosen octree level(s)
                if ~isempty(find(lvl == vis_levels_wavelet))
                    w_vis_flag = 1;
                    figure;
                    %Plot the input point cloud
                    scatter3(ptcloud(:, 1), ptcloud(:, 2), ptcloud(:, 3), 5, [ptcloud(:, 7)./255, ptcloud(:, 8)./255, ptcloud(:, 9)./255], 'filled');
                    axis equal; axis off;
                    hold on;
                    %Outline the current parent cell in thick black lines 
                    %(i.e., connect all the vertices in parent_corner_coords  
                    %to the other vertices in parent_corner_coords, with 
                    %which they share an edge). Since we are only plotting
                    %one cell, we know in advance how the vertices are 
                    %connected, so first construct a matrix of edge indices.
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
                    %Connect the vertices in the current parent cell, 
                    %according to edges_onecell
                    for pvtx1 = 1:8
                        for pvtx2 = 1:3
                            h_w(2) = plot3([parent_corner_coords(pvtx1, 1) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 1)], [parent_corner_coords(pvtx1, 2) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 2)], [parent_corner_coords(pvtx1, 3) parent_corner_coords(edges_onecell(pvtx1, pvtx2), 3)], 'k', 'LineWidth', 3);
                            hold on;
                        end
                    end
                    %Outline the current child cell in thick, dotted black 
                    %lines, by connecting the vertices in child_corner_coords 
                    %according to edges_onecell
                    for cvtx1 = 1:8
                        for cvtx2 = 1:3
                            h_w(3) = plot3([child_corner_coords(cvtx1, 1) child_corner_coords(edges_onecell(cvtx1, cvtx2), 1)], [child_corner_coords(cvtx1, 2) child_corner_coords(edges_onecell(cvtx1, cvtx2), 2)], [child_corner_coords(cvtx1, 3) child_corner_coords(edges_onecell(cvtx1, cvtx2), 3)], '--k', 'LineWidth', 3);
                            hold on;
                        end
                    end               
                    %Plot the current corner point (in red, with a black 
                    %outline) on the current figure
                    h_w(4) = scatter3(child_corner_coords(cnr, 1), child_corner_coords(cnr, 2), child_corner_coords(cnr, 3), 80, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
                    hold on;   
                end %End first part of the visualization code for the current corner   
                
                %---------------------------------------------------------%
                
                %For the current corner, check if it is on a parent edge 
                %(i.e., if it shares at least 2 same coordinates (out of x,
                %y, or z) with two of the parent vertices), or on a parent 
                %face (i.e., only 1 of its coordinates (either x, or y, or 
                %z) is the same as four parents' coordinates (must have the 
                %same coordinate in common in this case: either x, or y, or 
                %z))  
                parent_row_inds = [];
                [x_same, ~] = find(child_corner_coords(cnr, 1) == parent_corner_coords(:, 1));
                [y_same, ~] = find(child_corner_coords(cnr, 2) == parent_corner_coords(:, 2));
                [z_same, ~] = find(child_corner_coords(cnr, 3) == parent_corner_coords(:, 3));
                temp_cat = [x_same; y_same; z_same];
                if ((isempty(x_same) + isempty(y_same) + isempty(z_same)) == 2)
                    %Corner is on a parent face
                    on_p_face = 1;
                    %Get the row indices of the parent vertices on this 
                    %face
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
                if sum(ismember(parent_corner_coords, child_corner_coords(cnr, :), 'rows') > 0)
                    %Get the index of the current parent corner
                    %parent_cnr_index = find(ismember(parent_corner_coords, child_corner_coords(cnr, :), 'rows') > 0);
                    %Find the row index of this coordinate inside 
                    %corner_coords at the parent octree level
                    %parent_row_index = find(ismember(corner_coords{lvl}, parent_corner_coords(parent_cnr_index, :), 'rows') > 0);
                    %parent_row_index = find(ismember(corner_coords{lvl}, child_corner_coords(cnr, :), 'rows') > 0); 
                    parent_row_index = find(sum(abs(child_corner_coords(cnr, :) - corner_coords{lvl}), 2) == 0, 1); %Faster than "ismember"
                    %Find the control point index for the parent corner.
                    %Although the length of parent_row_index may sometimes
                    %be > 1, and so more than one control point index may
                    %be found below, in this case the result should still 
                    %be the same index, just repeated. So extract only the 
                    %unique control point index found below (should be just
                    %one).
                    %parent_ctrlpt_index = unique(ctrl_pts_pointers{lvl}(parent_row_index));
                    parent_ctrlpt_index = ctrl_pts_pointers{lvl}(parent_row_index);
                    %Do nothing (the signal on this vertex is a low-pass
                    %coefficient and has already been reconstructed),
                    %except insert a 0 here in the wavelet_coeffs cell 
                    %array, and transfer the reconstructed control point 
                    %over from the parent corner at the previous octree 
                    %level (we want the reconstructed_control_points values  
                    %and wavelet_coeffs values at each level to correspond
                    %to the same locations in unique_coords)
                    wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = 0;
                    reconstructed_control_points{lvl + 1}(possible_inds(cnr)) = reconstructed_control_points{lvl}(parent_ctrlpt_index);
                    if w_vis_flag == 1
                        hold off;
                        legend(h_w, 'Root cell (for reference)', ['Current parent cell (at level ' num2str(lvl) ')'], ['Current child cell (at level ' num2str(lvl + 1) ')'], 'Current child corner', 'Location', 'best');
                        title({'Computing Wavelet Coefficients', ['for Each Child of Each Occupied Octree Cell at Level ' num2str(lvl)]});
                        w_vis_flag = 0; %Reset flag for next octree cell
                    end
                    continue;
                %If this corner vertex lies on a parent edge
                elseif on_p_edge == 1
                    %Get the indices of the rows of these coordinates in
                    %parent_corner_coords (should be only 2 rows) 
                    if w_vis_flag == 1
                        for pc = 1:length(parent_row_inds)
                            %Circle the parent corners (in blue)
                            h_w(5) = scatter3(parent_corner_coords(parent_row_inds(pc), 1), parent_corner_coords(parent_row_inds(pc), 2), parent_corner_coords(parent_row_inds(pc), 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                            hold on;
                        end
                    end
                    %Get the Bezier control points stored at the corner 
                    %vertices defined by parent_row_inds
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt1 = reconstructed_control_points{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstructed_control_points{lvl}(ctrlpt2_ptr);
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
                    all_ctrlpts_ptrs = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpt1_ptr = all_ctrlpts_ptrs(parent_row_inds(1));
                    ctrlpt2_ptr = all_ctrlpts_ptrs(parent_row_inds(2));
                    ctrlpt3_ptr = all_ctrlpts_ptrs(parent_row_inds(3));
                    ctrlpt4_ptr = all_ctrlpts_ptrs(parent_row_inds(4));
                    ctrlpt1 = reconstructed_control_points{lvl}(ctrlpt1_ptr);
                    ctrlpt2 = reconstructed_control_points{lvl}(ctrlpt2_ptr);
                    ctrlpt3 = reconstructed_control_points{lvl}(ctrlpt3_ptr);
                    ctrlpt4 = reconstructed_control_points{lvl}(ctrlpt4_ptr);
                    %Average the signal (Bezier control points) on the 4 
                    %vertices of the parent face
                    avg_signal = (ctrlpt1 + ctrlpt2 + ctrlpt3 + ctrlpt4)/4;  
                %If this corner vertex lies in the centre of the parent
                %block (i.e., neither on a parent's edge nor on a parent's 
                %face)
                else                
                    if w_vis_flag == 1
                        for pc = 1:8
                            %Circle the parent corners (in blue)
                            h_w(5) = scatter3(parent_corner_coords(pc, 1), parent_corner_coords(pc, 2), parent_corner_coords(pc, 3), 80, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
                            hold on;
                        end
                    end
                    %Get the Bezier control points stored at each of the 8
                    %corner vertices of the parent cell
                    ctrlpts_pointers = ctrl_pts_pointers{lvl}((occ_cell*8 - 7):occ_cell*8, :);
                    ctrlpts = reconstructed_control_points{lvl}(ctrlpts_pointers);
                    %Average the signal (Bezier control points) on the 8
                    %vertices of the parent block
                    avg_signal = mean(ctrlpts);
                end
                %Get the Bezier control point stored at the current corner
                %of the current child cell
                child_ctrlpt = control_points{lvl + 1}(possible_inds(cnr));
                %Subtract the average signal from the signal (Bezier
                %control point) at the current child vertex. The result is 
                %the high-pass transform (wavelet) coefficient of the child 
                %vertex.
                wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = child_ctrlpt - avg_signal;
                %Quantize the wavelet coefficient computed above
                wavelet_coeffs{lvl + 1}(possible_inds(cnr)) = quantize_uniform_scalar(wavelet_coeffs{lvl + 1}(possible_inds(cnr)), q_stepsize);
                %Add the quantized wavelet coefficient to avg_signal, to 
                %obtain the reconstructed signal (control point) at the 
                %current child vertex 
                reconstructed_control_points{lvl + 1}(possible_inds(cnr)) = wavelet_coeffs{lvl + 1}(possible_inds(cnr)) + avg_signal;
                if w_vis_flag == 1
                    hold off;
                    legend(h_w, 'Root cell (for reference)', ['Current parent cell (at level ' num2str(lvl) ')'], ['Current child cell (at level ' num2str(lvl + 1) ')'], 'Current child corner', 'Parent corners', 'Location', 'best');
                    title({'Computing Wavelet Coefficients', ['for Each Child of Each Occupied Octree Cell at Level ' num2str(lvl)]});
                    w_vis_flag = 0; %Reset flag for next octree cell
                end   
            end %End corners
        end %End children
        %Increment parent_cnr_coords_cntr before moving on to a new 
        %occupied (parent) cell at the current octree level
        parent_cnr_coords_cntr = parent_cnr_coords_cntr + 8; 
    end %End occupied cells at current level
    %Arrange all the wavelet coefficients produced for the child octree 
    %level, into a column vector instead of a row vector (purely for 
    %visualization reasons: it is easier to scroll through a long column 
    %vector than a row vector). Also, wavelet_coeffs will then be in the 
    %same format as the control_points cell array.
    wavelet_coeffs{lvl + 1} = (wavelet_coeffs{lvl + 1})';
    %Also arrange all reconstructed control points at the next level, into
    %a column vector
    reconstructed_control_points{lvl + 1} = (reconstructed_control_points{lvl + 1})';
    wcfs_time = toc;
    disp(['Time taken to compute wavelet coefficients and reconstructed control points for level ' num2str(lvl + 1) ': ' num2str(wcfs_time) ' seconds']);
    disp('------------------------------------------------------------');
end %End current octree level
% waveletcfs_time = toc;
% disp(' ');
% disp('************************************************************');
% disp(['Time taken to compute all wavelet coefficients: ' num2str(waveletcfs_time) ' seconds']);
% disp('************************************************************');