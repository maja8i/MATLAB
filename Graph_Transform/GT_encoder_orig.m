%-------------------------------------------------------------------------%

%Encoder based on computing the Graph Transform for the geometry (point 
%locations) within each occupied cell of an octree, at a user-defined 
%octree level.

%Requires directory "Phil": add this directory plus its sub-directories to
%the current MATLAB path.

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

%p: Percentage of the largest spectral coefficients to use for point cloud 
%   reconstruction, from each occupied cell at octree level GT_block_lvl, 
%   for which the Graph Transform is computed. NOTE: p represents the
%   percentage of ALL spectral coefficients selected, not coefficients
%   corresponding to x, y, and z separately. p = 0 corresponds to no
%   coefficients being selected for reconstruction; p = 100 corresponds to
%   ALL coefficients being selected.

%q_stepsize: Step size used for uniform scalar quantization. Values for 
%            step_size should ideally be powers of 2. q_stepsize = 1 
%            corresponds to just rounding the input to the nearest integer, 
%            so the least amount of (uniform) quantization. q_stepsize = 0 
%            means that no quantization should be applied. 

%---- OUTPUTS ----

%nbr_occ_voxels_vec: The numbers of occupied voxels in each occupied cell 
%                    at level GT_block_lvl of the generated octree. 

%quantized_spectral_coeffs_cat: Quantized selected spectral coefficients, 
%                               arranged in a single column vector. 

%eig_indices_cat: Eigenvector indices corresponding to the selected
%                 spectral coefficients, also organized in a single column 
%                 vector.

%XYZ_indices_cat: Indices to indicate whether each of the selected spectral
%                 coefficients should be used to reconstruct the X, Y, or Z 
%                 input geometry vector (1 => X, 2 => Y, 3 => Z).

%lone_voxels_cat: Quantized (x, y, z) positions of lone voxels in octree  
%                 cells on which the Graph Transform could not be applied, 
%                 organized in a single column vector in the format: 
%                 [x1
%                  .
%                  .
%                  xM
%                  y1
%                  .
%                  .
%                  yM
%                  z1
%                  .
%                  .
%                  zM], where M is the total number of lone voxels across
%                  all occupied cells in the generated octree.

%rlgr_bits: Vector of bitrates for the transmitted data, computed using 
%           Adaptive Run-Length Golomb Rice encoding, arranged in a
%           2-column vector in the format: bpv_rlgr total_bits_rlgr 

%entropy_bits: Vector of bitrates for the transmitted data, computed using 
%              theoretical (practical) entropy, arranged in a 2-column 
%              vector in the format: bpv_entropy total_bits_entropy 

%-------------------------------------------------------------------------%

function [nbr_occ_voxels_vec, quantized_spectral_coeffs_cat, eig_indices_cat, XYZ_indices_cat, lone_voxels_cat, rlgr_bits, entropy_bits] = GT_encoder_orig(ptcloud_file, b, GT_block_lvl, p, q_stepsize)

disp('===================== ENCODER RUNNING ======================');

%------------------ Octree Construction and Processing -------------------%

%Construct octree of a given input point cloud
[myOT, mortonCodes_sorted, xyz_sorted] = construct_octree(ptcloud_file, b);

%Get the total number of occupied nodes (cells) at level GT_block_lvl of
%the octree myOT
total_occupied_cells = myOT.NodeCount(GT_block_lvl);

%Extract the set of occupied voxel coordinates (x, y, z) at all levels of
%the octree myOT - ****should modify this function so that it can compute
%voxel coordinates at any given level, without having to do all levels****
[~, occupied_voxel_coords] = extract_occupied_voxels(myOT, mortonCodes_sorted, xyz_sorted);
%Extract only the occupied voxel coordinates at level GT_block_lvl
occupied_voxel_coords_GT = cell(1, myOT.NodeCount(GT_block_lvl));
for i = 1:myOT.NodeCount(GT_block_lvl)
    occupied_voxel_coords_GT{1, i} = occupied_voxel_coords{GT_block_lvl, i};
end

%Store the numbers of occupied voxels in each occupied cell at level
%GT_block_lvl of myOT, in a separate array
nbr_occ_voxels_vec = myOT.DescendantCount{GT_block_lvl};

%---------------------------- Initializations ----------------------------%

%Initialize cell arrays to hold the adjacency matrix, degrees matrix, and
%Laplacian matrix for each occupied cell (node) of the octree myOT at level
%GT_block_lvl
A = cell(1, total_occupied_cells); %Adjacency matrices 
D = cell(1, total_occupied_cells); %Degrees matrices 
L = cell(1, total_occupied_cells); %Laplacian matrices 

%Initialize counter that keeps track of how many octree cells at
%GT_block_lvl have only 1 occupied voxel
one_occ_voxel_cnt = 0;

%Initialize a matrix to store the (x, y, z) positions of lone voxels in
%cells, on which the Graph Transform cannot be applied
lone_voxels = [];   %Will be an M x 3 matrix, where M is the total number of lone voxels at level GT_block_lvl
lone_voxels_quantized = []; %Quantized versions of the lone voxel positions
lone_voxel_cntr = 1;

%Initialize cell arrays to hold the eigenvector and eigenvalue matrices for
%each occupied cell (node) of the octree myOT at level GT_block_lvl
eigenvectors = cell(1, total_occupied_cells);
eigenvalues = cell(1, total_occupied_cells);

%Initialize a cell array to hold the spectral coefficients for each occupied
%cell (node) of the octree myOT at level GT_block_lvl
spectral_coeffs = cell(1, total_occupied_cells);

%Initialize a cell array to store the number of spectral coefficients
%selected in each occupied cell of the octree myOT at level GT_block_lvl
nbr_coeffs_selected = cell(1, total_occupied_cells);

%Initialize a cell array to store the indices of the ordered spectral
%coefficients in each occupied cell (node) of the octree myOT at level
%GT_block_lvl
I_coeffs_cell = cell(1, total_occupied_cells);

%Initialize a cell array to store the quantized set of selected spectral
%coefficients at each occupied cell (node) of the octree myOT at level
%GT_block_lvl
quantized_spectral_coeffs = cell(1, total_occupied_cells);

%------------------- Individual Octree Cell Processing -------------------%

%For each occupied cell at level GT_block_lvl ...
for n = 1:total_occupied_cells
    %Get the total number of occupied voxels in the current cell
    nbr_occ_voxels = nbr_occ_voxels_vec(n);
    
    %---------------------------- Lone Voxels ----------------------------%
    
    %If the current cell contains only 1 occupied voxel, we can't construct
    %a graph for this cell, so don't consider it for the Graph Transform, 
    %but record the corresponding voxel's (x, y, z) coordinates for
    %transmission to the decoder
    if nbr_occ_voxels == 1
        %Record this voxel's (x, y, z) coordinates
        lone_voxels(lone_voxel_cntr, :) = occupied_voxel_coords_GT{n};
        %Quantize the lone voxel positions for the current cell, using 
        %uniform scalar quantization
        lone_voxels_quantized(lone_voxel_cntr, 1) = quantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 1), q_stepsize);
        lone_voxels_quantized(lone_voxel_cntr, 2) = quantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 2), q_stepsize);
        lone_voxels_quantized(lone_voxel_cntr, 3) = quantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 3), q_stepsize);
        lone_voxel_cntr = lone_voxel_cntr + 1;
        %disp('------------------------------------------------------------');
        %disp(['Node ' num2str(n) ' contains only 1 voxel. Skipping ...']);
        one_occ_voxel_cnt = one_occ_voxel_cnt + 1;
        %disp(['one_occ_voxel_cnt = ' num2str(one_occ_voxel_cnt)]);
        continue;
    end
    
    %------------------------ A, D, L Construction -----------------------%

    %For a fully connected graph, vertex degrees will just be equal to the
    %number of occupied voxels in the cell minus 1
    vtx_degs = (nbr_occ_voxels - 1)*ones(nbr_occ_voxels, 1); 
    
    %Construct A, D, and L (combinatorial, unweighted) matrices
    A{n} = ~eye(nbr_occ_voxels);    %Adjacency matrix has zeros down the main diagonal and ones elsewhere
    D{n} = diag(vtx_degs, 0);   %Degrees matrix has vertex degrees down the main diagonal and zeros elsewhere
    L{n} = D{n} - A{n};   %Combinatorial Laplacian matrix has vertex degrees down the main diagonal and -1 elsewhere   
    
    %----------------- Eigenvector/Eigenvalue Computation ----------------%
    
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
        
    %Get the set of occupied voxel coordinates (x, y, z) in the current 
    %cell
    current_voxel_set = occupied_voxel_coords_GT{:, n}; %N x 3 matrix
    
    %--------------------- Graph Transform (Geometry) --------------------%
    
    %Project the geometry vectors (x, y, z) for the occupied voxels in the 
    %current cell, onto the matrix of sorted eigenvectors (basis)
    spectral_coeffs{n} = current_voxel_set'*eigvecs_sorted;
    
    %------------------- Spectral Coefficient Selection ------------------%
    
    %Select a certain percentage (p) of the largest-magnitude spectral 
    %coefficients to use for point cloud reconstruction. NOTE: p represents 
    %the percentage of ALL spectral coefficients selected for this cell, 
    %not coefficients corresponding to x, y, and z separately.
    [nrows, ncols] = size(spectral_coeffs{n});
    nbr_coeffs_selected{n} = ceil((p/100)*numel(spectral_coeffs{n}));
    [~, I_coeffs] = sort(abs(spectral_coeffs{n}(:)), 'descend');
    coeffs_sorted = spectral_coeffs{n}(I_coeffs);   %Keep original signs (+/-) of sorted coefficients
    coeffs_selected = coeffs_sorted(1:nbr_coeffs_selected{n});
    I_coeffs_cell{n} = I_coeffs;
    
    %----------------- Spectral Coefficient Quantization -----------------%
    
    %Quantize the selected set of spectral coefficients for the current
    %cell, using uniform scalar quantization
    quantized_spectral_coeffs{n} = quantize_uniform_scalar(coeffs_selected, q_stepsize);   
end
disp('------------------------------------------------------------');
disp(['Total no. of octree cells with 1 occupied voxel at level ' num2str(GT_block_lvl) ': ' num2str(one_occ_voxel_cnt)]);

%Arrange the occupied cells at level GT_block_lvl, in order of importance
%(most important cells first)? If used, the below will have to be arranged
%in the corresponding order and extra data will need to be encoded and
%transmitted to indicate the cell ordering.

%---------------- Distributions of Data to be Transmitted ----------------%

%No. of occupied voxels in each occupied cell at level GT_block_lvl (these
%are all unsigned, positive integers)
figure;
histogram(nbr_occ_voxels_vec);   
title({'Distribution of nbr\_occ\_voxels\_vec', '(No. of Occupied Voxels in Each Occupied Cell at the Chosen Octree Level)'});

%Quantized set of selected spectral coefficients (can be non-negative or
%negative integers)
quantized_spectral_coeffs_cat = [];
for n = 1:total_occupied_cells
    %Concatenate all selected coefficients, from all occupied octree cells,
    %into a single column vector
    quantized_spectral_coeffs_cat = [quantized_spectral_coeffs_cat; quantized_spectral_coeffs{n}];
end
figure;
histogram(quantized_spectral_coeffs_cat);   
if q_stepsize == 0
    title({'Distribution of quantized\_spectral\_coeffs\_cat', '(Unquantized Chosen Spectral Coefficients)'});
else 
    title({'Distribution of quantized\_spectral\_coeffs\_cat', ['(Chosen Spectral Coefficients, Quantized with Stepsize = ' num2str(q_stepsize) ')']});
end

%Eigenvector indices corresponding to selected spectral coefficients, 
%arranged in the corresponding order (all positive integers), and XYZ 
%indices that indicate which input geometry vector (X, Y, or Z) each
%transmitted spectral coefficient will be used to reconstruct (can have
%values of 1, 2, or 3)
eig_indices_cat = [];
XYZ_indices_cat = [];
for n = 1:total_occupied_cells
    %Compute the row and column indices of each selected spectral
    %coefficient in the current octree cell
    [chosen_coeff_rows, chosen_eig_inds] = ind2sub([3, size(spectral_coeffs{n}, 2)], I_coeffs_cell{n}(1:nbr_coeffs_selected{n}));
    %Concatenate all eigenvector indices corresponding to the selected 
    %spectral coefficients, from all occupied octree cells, into a single
    %column vector
    eig_indices_cat = [eig_indices_cat; chosen_eig_inds];
    %Concatenate all chosen_coeffs_rows, from all occupied octree cells,
    %into a single column vector. These indices indicate whether the  
    %corresponding spectral coefficient will be used to reconstruct the X,
    %Y, or Z vetor (1 => X, 2 => Y, 3 => Z). 
    XYZ_indices_cat = [XYZ_indices_cat; chosen_coeff_rows];
end
figure;
histogram(eig_indices_cat);   
title({'Distribution of eig\_indices\_cat', '(Eigenvector Indices Corresponding to Chosen Spectral Coefficients)'});
figure;
histogram(XYZ_indices_cat);   
title({'Distribution of XYZ\_indices\_cat', '(Indices Indicating which Geometry Vector (X, Y, or Z)', 'Each Spectral Coefficient will be Used to Reconstruct)'});

%Quantized lone voxels (can be non-negative or negative integers)
lone_voxels_cat = [lone_voxels_quantized(:, 1); lone_voxels_quantized(:, 2); lone_voxels_quantized(:, 3)];  %Concatenate columns to form a single-column vector
figure;
histogram(lone_voxels_cat);   
if q_stepsize == 0
    title({'Distribution of lone\_voxels\_cat', '(Unquantized Lone Voxel Coordinates)'});
else 
    title({'Distribution of lone\_voxels\_cat', ['(Lone Voxel Coordinates, Quantized with Stepsize = ' num2str(q_stepsize) ')']});
end

%-------------------- Encoding Data to be Transmitted --------------------%

%NOTES: rlgr() is Adaptive Run-Length Golomb Rice encoding; entropy()
%computes the theoretical (practical) Shannon entropy of the given sequence
%of symbols. 

%Encode the data to be transmitted to the decoder (i.e., figure out how
%many bits would be needed to encode this data) ...

%No. of occupied voxels in each occupied cell at level GT_block_lvl (these
%are all unsigned, positive integers)
bits_nbr_occ_voxels_rlgr = rlgr(nbr_occ_voxels_vec);    
disp('------------------------------------------------------------');
disp(['bits_nbr_occ_voxels_rlgr: ' num2str(bits_nbr_occ_voxels_rlgr)]);

bits_nbr_occ_voxels_entropy_persymbol = entropy_calc(nbr_occ_voxels_vec);   %Avg. no. of bits per symbol
bits_nbr_occ_voxels_entropy = bits_nbr_occ_voxels_entropy_persymbol*length(nbr_occ_voxels_vec); %Total no. of bits for all symbols
disp(['bits_nbr_occ_voxels_entropy: ' num2str(bits_nbr_occ_voxels_entropy) ' (' num2str(bits_nbr_occ_voxels_entropy_persymbol) ' bits per symbol)']);
disp(' ');

%Quantized set of selected spectral coefficients (can be non-negative or
%negative integers)
bits_quant_coeffs_rlgr = rlgr(quantized_spectral_coeffs_cat);
disp(['bits_quant_coeffs_rlgr: ' num2str(bits_quant_coeffs_rlgr)]);

bits_quant_coeffs_entropy_persymbol = entropy_calc(quantized_spectral_coeffs_cat);  %Avg. no. of bits per symbol
bits_quant_coeffs_entropy = bits_quant_coeffs_entropy_persymbol*length(quantized_spectral_coeffs_cat);  %Total no. of bits for all symbols
disp(['bits_quant_coeffs_entropy: ' num2str(bits_quant_coeffs_entropy) ' (' num2str(bits_quant_coeffs_entropy_persymbol) ' bits per symbol)']);
disp(' ');

%Eigenvector indices corresponding to selected spectral coefficients, 
%arranged in the corresponding order (all positive integers)
bits_eig_indices_rlgr = rlgr(eig_indices_cat);
disp(['bits_eig_indices_rlgr: ' num2str(bits_eig_indices_rlgr)]);

bits_eig_indices_entropy_persymbol = entropy_calc(eig_indices_cat); %Avg. no. of bits per symbol  
bits_eig_indices_entropy = bits_eig_indices_entropy_persymbol*length(eig_indices_cat);  %Total no. of bits for all symbols  
disp(['bits_eig_indices_entropy: ' num2str(bits_eig_indices_entropy) ' (' num2str(bits_eig_indices_entropy_persymbol) ' bits per symbol)']);
disp(' ');

%XYZ indices that indicate which input geometry vector (X, Y, or Z) each
%transmitted spectral coefficient will be used to reconstruct (can have
%values of 1, 2, or 3)
bits_XYZ_rlgr = rlgr(XYZ_indices_cat);
disp(['bits_XYZ_rlgr: ' num2str(bits_XYZ_rlgr)]);

bits_XYZ_entropy_persymbol = entropy_calc(XYZ_indices_cat); %Avg. no. of bits per symbol
bits_XYZ_entropy = bits_XYZ_entropy_persymbol*length(XYZ_indices_cat);  %Total no. of bits for all symbols  
disp(['bits_XYZ_entropy: ' num2str(bits_XYZ_entropy) ' (' num2str(bits_XYZ_entropy_persymbol) ' bits per symbol)']);
disp(' ');
    
%Quantized lone voxels (can be non-negative or negative integers)
bits_quant_lone_voxels_rlgr = rlgr(lone_voxels_cat);
disp(['bits_quant_lone_voxels_rlgr: ' num2str(bits_quant_lone_voxels_rlgr)]);

bits_quant_lone_voxels_entropy_persymbol = entropy_calc(lone_voxels_cat);   %Avg. no. of bits per symbol
bits_quant_lone_voxels_entropy = bits_quant_lone_voxels_entropy_persymbol*length(lone_voxels_cat);  %Total no. of bits for all symbols 
disp(['bits_quant_lone_voxels_entropy: ' num2str(bits_quant_lone_voxels_entropy) ' (' num2str(bits_quant_lone_voxels_entropy_persymbol) ' bits per symbol)']);
disp(' ');

%TOTAL BITS
total_bits_rlgr = bits_nbr_occ_voxels_rlgr + bits_quant_coeffs_rlgr + bits_eig_indices_rlgr + bits_XYZ_rlgr + bits_quant_lone_voxels_rlgr; 
disp('------------------------------------------------------------');
disp(['TOTAL BITS (RLGR): ' num2str(total_bits_rlgr)]);
total_bits_entropy = bits_nbr_occ_voxels_entropy + bits_quant_coeffs_entropy + bits_eig_indices_entropy + bits_XYZ_entropy + bits_quant_lone_voxels_entropy; 
disp(['TOTAL BITS (ENTROPY): ' num2str(total_bits_entropy)]);
disp(' ');

%Bits per voxel
bpv_rlgr = total_bits_rlgr/sum(nbr_occ_voxels_vec);
disp(['bpv_rlgr: ' num2str(bpv_rlgr)]);
bpv_entropy = total_bits_entropy/sum(nbr_occ_voxels_vec);
disp(['bpv_entropy: ' num2str(bpv_entropy)]);

%Store the bpv values and the total_bits values in their respective vectors
rlgr_bits = [bpv_rlgr, total_bits_rlgr];
entropy_bits = [bpv_entropy, total_bits_entropy];

disp('------------------------------------------------------------');
disp('Encoder finished.');
disp('============================================================');




