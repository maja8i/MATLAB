%-------------------------------------------------------------------------%

%Decoder corresponding to GT_encoder_orig.

%---- INPUTS ----

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

%q_stepsize: Step size used for uniform scalar quantization. Values for 
%            step_size should ideally be powers of 2. q_stepsize = 1 
%            corresponds to just rounding the input to the nearest integer, 
%            so the least amount of (uniform) quantization. q_stepsize = 0 
%            means that no quantization was applied at the encoder.

%p: Percentage of the largest spectral coefficients to use for point cloud 
%   reconstruction, from each occupied cell at octree level GT_block_lvl, 
%   for which the Graph Transform was computed. NOTE: p represents the
%   percentage of ALL spectral coefficients selected, not coefficients
%   corresponding to x, y, and z separately. p = 0 corresponds to no
%   coefficients being selected for reconstruction; p = 100 corresponds to
%   ALL coefficients being selected.

%---- OUTPUTS ----

%recon_xyz: Reconstructed geometry (x, y, z positions) of the input point 
%           cloud, arranged in an N x 3 matrix, where N is the total number 
%           of points in the input point cloud.   

%-------------------------------------------------------------------------%

function recon_xyz = GT_decoder_orig(nbr_occ_voxels_vec, quantized_spectral_coeffs_cat, eig_indices_cat, XYZ_indices_cat, lone_voxels_cat, q_stepsize, p)

disp('===================== DECODER RUNNING ======================');

%Work out the total number of occupied octree cells at the chosen level
total_occupied_cells = length(nbr_occ_voxels_vec);

%---------------------------- Initializations ----------------------------%

%Initialize cell arrays to hold the adjacency matrix, degrees matrix, and
%Laplacian matrix for each occupied cell (node) of the octree at the
%chosen level
A = cell(1, total_occupied_cells); %Adjacency matrices 
D = cell(1, total_occupied_cells); %Degrees matrices 
L = cell(1, total_occupied_cells); %Laplacian matrices

%Initialize cell arrays to hold the eigenvector and eigenvalue matrices for
%each occupied cell (node) of the octree at the chosen level
eigenvectors = cell(1, total_occupied_cells);
eigenvalues = cell(1, total_occupied_cells);

%Initialize a cell array to store the dequantized set of selected spectral
%coefficients at each occupied cell (node) of the octree at the chosen 
%level
dequantized_spectral_coeffs = cell(1, total_occupied_cells);

%Initialize cell arrays to hold the reconstructed voxel (x, y, z)
%positions for each occupied cell (node) at the chosen level of the octree
recon_xyz = cell(1, total_occupied_cells);

%------------------------------ Lone Voxels ------------------------------%

%Arrange all lone voxels into M x 3 format, where M is the total no. of
%lone voxels and each column represents the x, y, and z coordinates
%respectively
lone_voxels(:, 1) = lone_voxels_cat(1:(length(lone_voxels_cat)/3));
lone_voxels(:, 2) = lone_voxels_cat(((length(lone_voxels_cat)/3) + 1):(2*(length(lone_voxels_cat)/3)));
lone_voxels(:, 3) = lone_voxels_cat((2*(length(lone_voxels_cat)/3) + 1):length(lone_voxels_cat));

%Counter to keep track of the number of all lone voxels that have been 
%processed for reconstruction
lone_voxel_cntr = 1;
%Counter to keep track of the number of spectral coefficients that have
%been processed for reconstruction
coeffs_cntr = 1;

%------------------- Individual Octree Cell Processing -------------------%

%For each occupied cell at the chosen octree level ...
for n = 1:total_occupied_cells
    
    %---------------------------- Lone Voxels ----------------------------%
    
    %If the current cell contains only 1 occupied voxel, the reconstructed
    %(x, y, z) locations for this voxel are just obtained from lone_voxels
    if nbr_occ_voxels_vec(n) == 1
        %Dequantize the lone voxel positions for the current cell
        recon_xyz{n}(1, 1) = dequantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 1), q_stepsize);    %x
        recon_xyz{n}(1, 2) = dequantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 2), q_stepsize);    %y
        recon_xyz{n}(1, 3) = dequantize_uniform_scalar(lone_voxels(lone_voxel_cntr, 3), q_stepsize);    %z
        lone_voxel_cntr = lone_voxel_cntr + 1;
        continue;
    end
    
    %------------------------ A, D, L Construction -----------------------%
    
    %For a fully connected graph, vertex degrees will just be equal to the
    %number of occupied voxels in the cell minus 1
    vtx_degs = (nbr_occ_voxels_vec(n) - 1)*ones(nbr_occ_voxels_vec(n), 1);
    
    %Construct A, D, and L (combinatorial, unweighted) matrices
    A{n} = ~eye(nbr_occ_voxels_vec(n));    %Adjacency matrix has zeros down the main diagonal and ones elsewhere
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
        
    %---------------------- Reconstruction (Geometry) --------------------%
    
    %Work out how many spectral coefficients were selected for the current
    %octree cell
    nbr_coeffs_selected = ceil((p/100)*(3*length(eigenvalues{n})));
    
    %Extract the spectral coefficients corresponding to the current octree
    %cell
    current_spectral_coeffs = quantized_spectral_coeffs_cat(coeffs_cntr:(coeffs_cntr + nbr_coeffs_selected - 1));
    
    %Dequantize the spectral coefficients corresponding to the current
    %octree cell
    dequant_current_spectral_coeffs = dequantize_uniform_scalar(current_spectral_coeffs, q_stepsize);
    
    %Extract the eigenvector indices corresdponding to the current octree
    %cell
    current_eig_inds = eig_indices_cat(coeffs_cntr:(coeffs_cntr + nbr_coeffs_selected - 1));
    
    %Extract the X/Y/Z indices corresponding to the current octree cell
    current_xyz_inds = XYZ_indices_cat(coeffs_cntr:(coeffs_cntr + nbr_coeffs_selected - 1));
        
    %Initialize the current location of the dequantized_spectral_coeffs
    %cell array to contain all zeros
    dequantized_spectral_coeffs{n} = zeros(3, length(eigenvalues{n}));
    
    %Insert the chosen spectral coefficient values into their
    %corresponding locations inside dequantized_spectral_coeffs{n}, based
    %on their eigenvector index values and XYZ index values
    for i = 1:length(current_xyz_inds)
        dequantized_spectral_coeffs{n}(current_xyz_inds(i), current_eig_inds(i)) = dequant_current_spectral_coeffs(i);
    end
    
    %Reconstruct the x, y, z locations of the occupied voxels in the
    %current octree cell
    recon_xyz_temp = dequantized_spectral_coeffs{n}*eigenvectors{n}';   
    recon_xyz{n} = recon_xyz_temp'; %N x 3 matrix
    
    %Update the counter, ready for processing the next octree cell
    coeffs_cntr = coeffs_cntr + nbr_coeffs_selected;  
end

disp('------------------------------------------------------------');
disp('Decoder finished.');
disp('============================================================');

