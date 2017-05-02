%-------------------------------------------------------------------------%

%Compute the Graph Transform for the geometry (point locations) within each
%occupied cell of an octree, at a user-defined octree level.

%Requires directories "Phil", "Quantization": add these directories plus
%their sub-directories to the current MATLAB path.

%---- INPUTS ----

%ptcloud_file: Name of the input PLY file, including the full path to this
%              file if it is not in the current MATLAB directory, e.g., 
%              ptcloud_file =
%              '\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10\boxer_voxelized10.ply'.

%b: Bit depth for Morton codes and octree. b also determines the number of
%   levels in the octree that will be generated (apart from the root
%   level). The total number of octree levels INCLUDING the root level will
%   therefore be: b + 1. IMPORTANT: If using a voxelized point cloud as 
%   input, b must be equal to the voxelization level used to produce that 
%   point cloud, e.g., b = 10 for voxelized10 point clouds, b = 11 for 
%   voxelized11 clouds, etc.

%GT_block_lvl: Desired octree level, at which the Graph Transform will be 
%              computed for each occupied cell. IMPORTANT: GT_block_lvl 
%              must be <= b.  

%L_type: Type of Laplacian to use for the Graph Transform computation.
%        Currently only considering combinatorial ("comb") Laplacians. 
%
%        Choices are:
%
%        'comb_unweighted': L = D - A, where A contains a weight of 1 for 
%                         each edge and zeros down the main diagonal; D has 
%                         vertex degrees down the main diagonal and zeros 
%                         elsewhere; so L := deg(vi) if i = j; -1 if i ~= j 
%                         and vi is adjacent to vj; 0 otherwise.
%
%        'comb_weighted': L = D - A, where A has a weight of 1/Euclidean 
%                       distance for each edge and zeros down the main 
%                       diagonal; D contains the sum of weights (rows or 
%                       columns in A) down the main diagonal and zeros 
%                       elsewhere; so L := sum(weights) if i = j; 
%                       -1/Euclidean distance if i ~= j and vi is adjacent
%                       to vj; 0 otherwise. Euclidean distance is the 
%                       Euclidean distance between pairs of points that are
%                       considered to be connected by an edge in the graph
%                       that A describes.

%p: Percentage of the largest spectral coefficients to use for point cloud 
%   reconstruction, from each occupied cell at octree level GT_block_lvl 
%   for which the Graph Transform is computed. NOTE: p represents the
%   percentage of ALL spectral coefficients selected, not coefficients
%   corresponding to x, y, and z separately. p = 0 corresponds to no
%   coefficients being selected for reconstruction; p = 100 corresponds to
%   ALL coefficients being selected.

%q_stepsize: Step size used for uniform scalar quantization. Values for 
%            step_size should ideally be powers of 2. q_stepsize = 1 
%            corresponds to just rounding the input to the nearest integer, 
%            so the least amount of (uniform) quantization.  

%---- OUTPUTS ----

%myOT: Octree class, with a number of different properties.

%occupied_voxel_coords_GT: Locations (x, y, z) of the occupied voxels in 
%                          all occupied cells at level GT_block_lvl of the 
%                          octre myOT.

%A: Adjacency matrices for all occupied cells (nodes) of the octree myOT at
%   level GT_block_lvl.

%D: Degrees matrices for all occupied cells (nodes) of the octree myOT at
%   level GT_block_lvl.

%L: Laplacian matrices for all occupied cells (nodes) of the octree myOT at
%   level GT_block_lvl.

%eigenvectors: Matrices of eigenvectors for each occupied cell (node) of
%              the octree myOT at level GT_block_level, organised in order 
%              of increasing eigenvalue magnitude.

%eigenvalues: Vectors of eigenvalues for each occupied cell (node) of the
%             octree myOT at level GT_block_level, organised in order 
%             of increasing magnitude.

%spectral_coeffs: Matrices of spectral coefficients for each occupied cell
%                 (node) of the octree myOT at level GT_block_lvl, 
%                 resulting from projecting the geometry vectors (x, y, z) 
%                 of the occupied voxels in each cell onto the matrix of 
%                 sorted eigenvectors.

%spectral_coeffs_thresholded: Thresholded version of spectral_coeffs, where
%                             only the selected coefficients maintain their
%                             original values, while all other spectral 
%                             coefficients are set to 0.

%quantized_spectral_coeffs: Quantized versions of the thresholded spectral
%                           coefficients.

%dequantized_spectral_coeffs: Reconstructed (dequantized) versions of the
%                             thresholded set of spectral coefficients.

%recon_xyz: Reconstructed voxel x, y, z positions for each occupied cell 
%           (node) at level GT_block_lvl of the octree myOT, when no
%           quantization has been used for the selected spectral 
%           coefficients. 

%recon_error: Difference (just subtraction) between each reconstructed 
%             (x, y, z) and the corresponding input voxel (x, y, z), for 
%             each occupied cell at level GT_block_lvl of the octree myOT,
%             when no quantization has been used for the selected spectral 
%             coefficients. 

%recon_xyz_quant: Voxel (x, y, z) positions for each occupied cell (node)
%                 at level GT_block_lvl of the octree myOT, reconstructed 
%                 from the dequantized set of selected spectral 
%                 coefficients.

%recon_quant_error: Difference (just subtraction) between each 
%                   reconstructed (x, y, z) and the corresponding input 
%                   voxel (x, y, z), for each occupied cell at level 
%                   GT_block_lvl of the octree myOT, when the reconstructed
%                   positions have been obtained by using a dequantized set
%                   of selected spectral coefficients.

%-------------------------------------------------------------------------%

function [myOT, occupied_voxel_coords_GT, A, D, L, eigenvectors, eigenvalues, spectral_coeffs, spectral_coeffs_thresholded, quantized_spectral_coeffs, dequantized_spectral_coeffs, recon_xyz, recon_error, recon_xyz_quant, recon_quant_error] = graph_transform_geom(ptcloud_file, b, GT_block_lvl, L_type, p, q_stepsize)

%Construct octree
[myOT, mortonCodes_sorted, xyz_sorted] = construct_octree(ptcloud_file, b);

%Extract the set of occupied voxel coordinates (x, y, z) at all levels of
%the octree myOT
[~, occupied_voxel_coords] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted);
%Extract only the occupied voxel coordinates at level GT_block_lvl
occupied_voxel_coords_GT = cell(1, myOT.NodeCount(GT_block_lvl));
for i = 1:myOT.NodeCount(GT_block_lvl)
    occupied_voxel_coords_GT{1, i} = occupied_voxel_coords{GT_block_lvl, i};
end

%Initialize cell arrays to hold the adjacency matrix, degrees matrix, and
%Laplacian matrix for each occupied cell (node) of the octree myOT at level
%GT_block_lvl
A = cell(1, myOT.NodeCount(GT_block_lvl)); %Adjacency matrices 
D = cell(1, myOT.NodeCount(GT_block_lvl)); %Degrees matrices 
L = cell(1, myOT.NodeCount(GT_block_lvl)); %Laplacian matrices 

%Initialize cell arrays to hold the eigenvector and eigenvalue matrices for
%each occupied cell (node) of the octree myOT at level GT_block_lvl
eigenvectors = cell(1, myOT.NodeCount(GT_block_lvl));
eigenvalues = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize cell arrays to hold the spectral coefficients for each occupied
%cell (node) of the octree myOT at level GT_block_lvl. spectral_coeffs will
%store all of the spectral coefficients that are produced from the Graph
%Transform in each cell; spectral_coeffs_thresholded will contain actual
%values for only the spectral coefficients that will be used in the point 
%cloud reconstruction, while the other coefficients will have values of 0.
spectral_coeffs = cell(1, myOT.NodeCount(GT_block_lvl));
spectral_coeffs_thresholded = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize a cell array to store the quantized set of selected spectral
%coefficients at each occupied cell (node) of the octree myOT at level
%GT_block_lvl
quantized_spectral_coeffs = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize a cell array to store the dequantized set of selected spectral
%coefficients at each occupied cell (node) of the octree myOT at level
%GT_block_lvl
dequantized_spectral_coeffs = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize cell arrays to hold the reconstructed voxel (x, y, z)
%positions for each occupied cell (node) at level GT_block_lvl of the
%octree myOT. recon_xyz will contain the reconstructed positions when no
%quantization has been used for the selected spectral coefficients; 
%recon_xyz_quant will contain the reconstructed positions when quantization
%has been used.
recon_xyz = cell(1, myOT.NodeCount(GT_block_lvl));
recon_xyz_quant = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize cell arrays to hold the reconstruction errors (obtained from
%straight subtraction) between each reconstructed (x, y, z) and the
%corresponding input voxel (x, y, z), for each occupied cell at level 
%GT_block_lvl of the octree myOT. recon_error will contain the errors for
%the case when no quantization has been used for the selected spectral
%coefficients; recon_quant_error will contain the errors for the case when
%the reconstructed positions have been obtained by using a dequantized set
%of selected spectral coefficients.  
recon_error = cell(1, myOT.NodeCount(GT_block_lvl));
recon_quant_error = cell(1, myOT.NodeCount(GT_block_lvl));

%Initialize counter that keeps track of how many octree cells at
%GT_block_lvl have only 1 occupied voxel
one_occ_voxel_cnt = 0;

%For each occupied octree cell at level GT_block_lvl, construct a fully
%connected graph of the occupied voxels in that cell (i.e., where every
%vertex is connected to every other vertex). Then, apply the Graph
%Transform in that cell, by using the eigenvectors of the Laplacian matrix 
%L_type of the graph in that cell as the basis functions.
for n = 1:myOT.NodeCount(GT_block_lvl)
    
    %Get the total number of occupied voxels in the current cell
    nbr_occ_voxels = myOT.DescendantCount{GT_block_lvl}(n);
    
    %Some occupied octree cells may contain only 1 voxel. We can't
    %construct a graph for these cells, so don't consider them. 
    if nbr_occ_voxels == 1
        disp('------------------------------------------------------------');
        disp(['Node ' num2str(n) ' contains only 1 voxel. Skipping ...']);
        one_occ_voxel_cnt = one_occ_voxel_cnt + 1;
        disp(['one_occ_voxel_cnt = ' num2str(one_occ_voxel_cnt)]);
        continue;
    end
            
    %For a fully connected graph, vertex degrees will just be equal to the
    %number of occupied voxels in the cell minus 1
    vtx_degs = (nbr_occ_voxels - 1)*ones(nbr_occ_voxels, 1);   
    
    %Get the set of occupied voxel coordinates (x, y, z) in the current 
    %cell
    current_voxel_set = occupied_voxel_coords_GT{:, n}; %N x 3 matrix
    
    %Construct A, D, and L matrices
    switch L_type
        case 'comb_unweighted'
            A{n} = ~eye(nbr_occ_voxels);    %Adjacency matrix has zeros down the main diagonal and ones elsewhere
            D{n} = diag(vtx_degs, 0);   %Degrees matrix has vertex degrees down the main diagonal and zeros elsewhere
            L{n} = D{n} - A{n};   %Combinatorial Laplacian matrix has vertex degrees down the main diagonal and -1 elsewhere   
        case 'comb_weighted'
            %Construct weighted graph, using weight = 1/(Euclidean distance 
            %between 'connected' points)
            A_temp = zeros(nbr_occ_voxels);    
            for i = 1:nbr_occ_voxels
                %Compute weights for lower triangular part of A_temp only,
                %since the matrix will be symmetric
                for j = 1:i
                    if i == j
                        continue;
                    else
                        euclid_dist = sqrt((current_voxel_set(i, 1) - current_voxel_set(j, 1))^2 + (current_voxel_set(i, 2) - current_voxel_set(j, 2))^2 + (current_voxel_set(i, 3) - current_voxel_set(j, 3))^2);
                        A_temp(i, j) = 1/euclid_dist;
                        %Copy the same weight value to the location (j, i),
                        %since A_temp is a symmetric matrix
                        A_temp(j, i) = A_temp(i, j);
                    end
                end
            end
            A{n} = A_temp;  %Adjacency matrix has zeros down the main diagonal and 1/Euclidean distance elsewhere
            D{n} = diag(sum(A{n}, 1), 0);  %Degrees matrix has sum of weights (rows or columns in adjacency matrix) down the main diagonal and zeros elsewhere
            L{n} = D{n} - A{n};   %Combinatorial Laplacian matrix has sum of weights down the main diagonal and -1/Euclidean distance elsewhere
    end
    
    %Compute all the eigenvectors and eigenvalues of L
    [eigvecs, eigvals_matrix] = eig(L{n});  %eigvecs are unit norm
    %eigvals_matrix is a square matrix with eigenvalues on the main 
    %diagonal: extract these into a column vector
    eigvals = diag(eigvals_matrix);

    %Sort the eigenvalues in order of increasing magnitude (using absolute
    %values, although in this case, since L is real and symmetric, the
    %eigenvalues should all be non-negative)
    [~, I_eigvals] = sort(abs(eigvals));
    eigvals_sorted = eigvals(I_eigvals);    %Keep original signs (+/-) of sorted eigenvalues
    %Sort the corresponding eigenvectors in the same order
    eigvecs_sorted = eigvecs(:, I_eigvals);

    %Store the sorted eigenvectors and eigenvalues in their respective cell
    %arrays
    eigenvectors{n} = eigvecs_sorted;
    eigenvalues{n} = eigvals_sorted;

    %Project the geometry vectors (x, y, z) for the occupied voxels in the 
    %current cell, onto the matrix of sorted eigenvectors (basis)
    spectral_coeffs{n} = current_voxel_set'*eigvecs_sorted;   
end

disp('------------------------------------------------------------');
disp(['Total no. of nodes with 1 occupied voxel: ' num2str(one_occ_voxel_cnt)]);

%In each occupied cell at level GT_block_lvl of myOT, select a certain 
%percentage (p) of spectral coefficients to transmit to the decoder, 
%quantize these coefficients, and then reconstruct the corresponding input
%voxel x, y, z positions.
for n = 1:myOT.NodeCount(GT_block_lvl)
    %If the Graph Transform was computed for the current cell (i.e., if the
    %Laplacian matrix is not empty due to there being only 1 point in this
    %cell) ...
    if ~isempty(L{n})
        %Select a certain percentage (p) of the largest-magnitude spectral 
        %coefficients to use for point cloud reconstruction. NOTE: p 
        %represents the percentage of ALL spectral coefficients selected 
        %for this cell, not coefficients corresponding to x, y, and z 
        %separately. 
        [nrows, ncols] = size(spectral_coeffs{n});
        nbr_coeffs_selected = ceil((p/100)*numel(spectral_coeffs{n}));
        [~, I_coeffs] = sort(abs(spectral_coeffs{n}(:)), 'descend');
        coeffs_sorted = spectral_coeffs{n}(I_coeffs);   %Keep original signs (+/-) of sorted coefficients
        %Initialize all locations of spectral_coeffs_thresholded to contain
        %zeros
        spectral_coeffs_thresholded{n} = zeros(nrows, ncols);
        %Only the selected spectral coefficients will maintain their 
        %original values in spectral_coeffs_thresholded; all other spectral 
        %coefficients will have been set to 0, above.
        coeffs_selected = coeffs_sorted(1:nbr_coeffs_selected);
        spectral_coeffs_thresholded{n}(I_coeffs(1:nbr_coeffs_selected)) = coeffs_selected;

        %Reconstruct the point (x, y, z) locations in the current cell, by 
        %by using the thresholded set of spectral coefficients
        recon_xyz_temp = spectral_coeffs_thresholded{n}*eigenvectors{n}'; 
        recon_xyz{n} = recon_xyz_temp'; %N x 3 matrix

        %Compute the reconstruction error (just by subtraction) between the
        %reconstructed voxel positions obtained above, and their 
        %corresponding original positions
        recon_error{n} = recon_xyz{n} - occupied_voxel_coords_GT{n};

        %Quantize the selected set of spectral coefficients for the current
        %cell, using uniform scalar quantization for x, y, and z components
        %separately
        quantized_spectral_coeffs{n}(1, :) = quantize_uniform_scalar(spectral_coeffs_thresholded{n}(1, :), q_stepsize);
        quantized_spectral_coeffs{n}(2, :) = quantize_uniform_scalar(spectral_coeffs_thresholded{n}(2, :), q_stepsize);
        quantized_spectral_coeffs{n}(3, :) = quantize_uniform_scalar(spectral_coeffs_thresholded{n}(3, :), q_stepsize);

        %Dequantize (reconstruct) the selected set of spectral coefficients 
        %from their corresponding quantized versions
        dequantized_spectral_coeffs{n}(1, :) = dequantize_uniform_scalar(quantized_spectral_coeffs{n}(1, :), q_stepsize);
        dequantized_spectral_coeffs{n}(2, :) = dequantize_uniform_scalar(quantized_spectral_coeffs{n}(2, :), q_stepsize);
        dequantized_spectral_coeffs{n}(3, :) = dequantize_uniform_scalar(quantized_spectral_coeffs{n}(3, :), q_stepsize);

        %Reconstruct the point (x, y, z) locations in the current octree 
        %cell, by using the dequantized set of spectral coefficients 
        recon_xyz_quant_temp = dequantized_spectral_coeffs{n}*eigenvectors{n}';   
        recon_xyz_quant{n} = recon_xyz_quant_temp'; %N x 3 matrix

        %Compute the reconstruction error (just by subtraction) for the
        %current cell, between the reconstructed voxel positions resulting
        %from using the dequantized set of selected spectral coefficients, 
        %and their corresponding original positions
        recon_quant_error{n} = recon_xyz_quant{n} - occupied_voxel_coords_GT{n};
    else
        %If the Graph Transform was not computed for the current cell
        %(because the cell contains only one point), just quantize the 
        %input voxel x, y, and z positions in these cells 
        current_voxel_set = occupied_voxel_coords_GT{:, n}; %N x 3 matrix
        quantized_voxels_x = quantize_uniform_scalar(current_voxel_set(:, 1), q_stepsize);
        quantized_voxels_y = quantize_uniform_scalar(current_voxel_set(:, 2), q_stepsize);
        quantized_voxels_z = quantize_uniform_scalar(current_voxel_set(:, 3), q_stepsize);
        
        %Reconstruct (dequantize) the quantized voxel positions
        dequantized_voxels_x = dequantize_uniform_scalar(quantized_voxels_x, q_stepsize);
        dequantized_voxels_y = dequantize_uniform_scalar(quantized_voxels_y, q_stepsize);
        dequantized_voxels_z = dequantize_uniform_scalar(quantized_voxels_z, q_stepsize);
        
        %Pass the dequantized voxel positions to the reconstruction matrix
        recon_xyz_quant{n}(1, :) = dequantized_voxels_x;
        recon_xyz_quant{n}(2, :) = dequantized_voxels_y;
        recon_xyz_quant{n}(3, :) = dequantized_voxels_z;
        
        %Compute the reconstruction error (just by subtraction) for the
        %current cell, between the reconstructed voxel positions and their
        %corresponding original positions
        recon_quant_error{n} = recon_xyz_quant{n} - occupied_voxel_coords_GT{n};
    end
end

    


