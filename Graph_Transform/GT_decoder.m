%-------------------------------------------------------------------------%

%Decoder corresponding to GT_encoder.

%---- INPUTS ----

%nbr_occ_voxels_diff_vec: Differences between the numbers of occupied
%                         voxels in successive occupied cells at level 
%                         GT_block_lvl of the generated octree. 

%thresh_spectral_coeffs_sorted: All quantized AC spectral coefficients, 
%                               after thresholding and sorting by 
%                               eigenvalue (smallest -> largest eigenvalue 
%                               magnitude across ALL occupied octree cells 
%                               at level GT_block_lvl). This vector is of
%                               size 3M x 1 and organized in the format:
%                               [x1 y1 z1 ... xM yM zM]', where M is the  
%                               total number of non-DC eigenvectors.

%DC_spectral_coeffs: Quantized DC spectral coefficients (those correspon-
%                    ding to eigenvalue of 0). The vector 
%                    DC_spectral_coeffs is of size 3D x 1, where D is the 
%                    total no. of 0 eigenvalues across all octree cells. It
%                    is arranged in the format: [x1 y1 z1 ... xD yD zD]'.

%q_stepsize: Step size used for uniform scalar quantization. Values for 
%            step_size should ideally be powers of 2. q_stepsize = 1 
%            corresponds to just rounding the input to the nearest integer, 
%            so the least amount of (uniform) quantization. q_stepsize = 0 
%            means that no quantization should be applied.

%A_construction: Option that tells the program how to construct the
%                adjacency matrix. Current choices are:
%
%                'fully_connected': every vertex in the octree cell is
%                connected to every other vertex in the same cell
%
%                'triangle_soup': the topology is defined by the triangles
%                in a triangle soup constructed at a chosen octree level

%---- OUTPUTS ----

%recon_xyz: Reconstructed geometry (x, y, z positions) of the input point 
%           cloud, arranged in an N x 3 matrix, where N is the total number 
%           of points in the input point cloud.   

%-------------------------------------------------------------------------%

function recon_xyz = GT_decoder(nbr_occ_voxels_diff_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, q_stepsize)

disp('===================== DECODER RUNNING ======================');

%Convert the difference vector nbr_occ_voxels_diff_vec to a vector
%containing the actual numbers of occupied voxels in each occupied cell at
%the chosen octree level
nbr_occ_voxels_vec = zeros(numel(nbr_occ_voxels_diff_vec), 1);
for i = 1:numel(nbr_occ_voxels_diff_vec)
    if i == 1
        nbr_occ_voxels_vec(i) = nbr_occ_voxels_diff_vec(1);
    else
        nbr_occ_voxels_vec(i) = nbr_occ_voxels_diff_vec(i) + nbr_occ_voxels_vec(i - 1);
    end
end

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

%Initialize a cell array to hold all the spectral coefficients (DC and AC),
%for each cell of the octree, at the chosen octree level, after these
%coefficients have been dequantized
spectral_coeffs = cell(1, total_occupied_cells);

%Initialize a cell array to hold the reconstructed x, y, z positions of the
%occupied voxels in each cell of the octree, at the chosen level
recon_xyz_cell = cell(1, total_occupied_cells);

%------------------- Individual Octree Cell Processing -------------------%

%For each occupied cell at the chosen octree level ...
for n = 1:total_occupied_cells
    
    %------------------------ A, D, L Construction -----------------------%
    
    %For a fully connected graph, vertex degrees will just be equal to the
    %number of occupied voxels in the cell minus 1
    vtx_degs = (nbr_occ_voxels_vec(n) - 1)*ones(nbr_occ_voxels_vec(n), 1);
    
    %Construct A, D, and L (combinatorial, unweighted) matrices
    switch A_construction
        case 'fully_connected'
            A{n} = ~eye(nbr_occ_voxels);    %Adjacency matrix has zeros down the main diagonal and ones elsewhere
            D{n} = diag(vtx_degs, 0);   %Degrees matrix has vertex degrees down the main diagonal and zeros elsewhere
        case 'triangle_soup'
            %Get the indices of the coarse-level faces that belong to the
            %current occupied octree cell
            face_inds = find(faceblocks == unique_faceblocks(n));
            %Extract the faces at indices face_inds
            current_faces = faces(face_inds);
            %Get the edges corresponding to current_faces
            E = [current_faces([1 2], :) current_faces([2 3], :) current_faces([3 1], :)];
            %Construct the adjacency matrix corresponding to edges E
            A = sparse(E(1, :), E(2, :), ones(size(E, 2), 1));
            A{n} = max(A, A');
            %Construct the degrees matrix corresponding to A{n}. The matrix
            %D{n} will have vertex degrees down the main diagonal and zeros
            %elsewhere.
            D{n} = spdiags(sum(A{n})', 0, nbr_occ_voxels, nbr_occ_voxels);
            %Expand A{n} and D{n} to full (non-sparse) matrices
            A{n} = full(A{n});
            D{n} = full(D{n});
    end
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
end  

%------------------------------- Dequantizing ----------------------------%

%Dequantize each of the values inside vector thresh_spectral_coeffs_sorted
dequant_thresh_spectral_coeffs_sorted = dequantize_uniform_scalar(thresh_spectral_coeffs_sorted, q_stepsize);

%Arrange the 3M x 1 vector dequant_thresh_spectral_coeffs_sorted into a 
%3 x M matrix
dequant_thresh_spectral_coeffs_sorted = reshape(dequant_thresh_spectral_coeffs_sorted, [3, length(dequant_thresh_spectral_coeffs_sorted)/3]);

%Dequantize each of the values inside vector DC_spectral_coeffs
dequant_DC_spectral_coeffs = dequantize_uniform_scalar(DC_spectral_coeffs, q_stepsize);

%Arrange the 3D x 1 vector dequant_DC_spectral_coeffs into a 3 x D matrix
dequant_DC_spectral_coeffs = reshape(dequant_DC_spectral_coeffs, [3, length(dequant_DC_spectral_coeffs)/3]);

%---------------------------- Eigenvalue Sorting -------------------------%

%Sort the eigenvalues computed for ALL occupied octree cells at level
%GT_block_lvl, in order from smallest to largest magnitude
all_eigenvalues = cell2mat(eigenvalues'); %N x 1 vector
%Find all "very small" eigenvalues and set them to 0. This is because these
%eigenvalues are effectively 0 anyway, but sometimes MATLAB produces very
%small negative or positive values that should be 0. Once they've all been
%set to 0, sorting them should just produce DC coefficients in the same
%order in which the octree cells were processed.
I_0eig = find(all_eigenvalues < 0.00000001);
all_eigenvalues(I_0eig) = 0;
[~, I_eigenvalues] = sort(abs(all_eigenvalues));

%-------------------------- Unsorting Received Data ----------------------%

% %Unnormalize each DC spectral coefficient and put it in its correct place
% %in spectral_coeffs, depending on which octree cell it belongs to 
% for n = 1:numel(spectral_coeffs)
%     spectral_coeffs{n}(:, 1) = dequant_DC_spectral_coeffs(:, n).*sqrt(nbr_occ_voxels_vec(n));   
% end

%Put each unquantized DC spectral coefficient in its correct place in
%spectral_coeffs, depending on which octree cell it belongs to
for n = 1:numel(spectral_coeffs)
    spectral_coeffs{n}(:, 1) = dequant_DC_spectral_coeffs(:, n);   
end

%Figure out which octree cell each of the sorted AC eigenvalues (and thus
%the corresponding column in dequant_thresh_spectral_coeffs_sorted) belongs
%to
AC_eigenvalues_inds = I_eigenvalues((total_occupied_cells + 1):end);
%Initialize an array to store the indices of the octree cells, to which the
%sorted AC eigenvalues belong
ACeig_cell_inds = zeros(1, length(AC_eigenvalues_inds));
for i = 1:length(AC_eigenvalues_inds)
    %Initialize counter that will keep track of how many occupied voxels
    %have been considered so far
    occ_vox_cntr = 0;
    %Get the unsorted eigenvalue index for the current eigenvalue in the
    %sorted set
    curr_unsort_eig_ind = AC_eigenvalues_inds(i);
    %Figure out which octree cell this eigenvalue belongs to
    for j = 1:length(nbr_occ_voxels_vec)
        %Number of occupied voxels considered so far
        occ_vox_cntr = occ_vox_cntr + nbr_occ_voxels_vec(j);
        if (curr_unsort_eig_ind <= occ_vox_cntr)
            %The eigenvalue belongs to the current octree cell, so record
            %this cell's index
            ACeig_cell_inds(i) = j;
            break;
        end
    end    
end

%Place all the AC spectral coefficients in their correct place in 
%spectral_coeffs, depending on which octree cell they belong to
for n = 1:numel(spectral_coeffs)
    %Extract the sorted eigenvalue indices, which belong to the current 
    %octree cell
    current_AC_eig_inds = find(ACeig_cell_inds == n);
    %Extract the corresponding AC spectral coefficients, from the sorted
    %set, and put them in their correct place in the spectral_coeffs cell 
    %array. If there are no AC coefficients (i.e., if the corresponding
    %octree cell has only 1 occupied voxel), then skip this step.
    if ~isempty(current_AC_eig_inds)
        spectral_coeffs{n}(:, 2:(length(current_AC_eig_inds)+1)) = dequant_thresh_spectral_coeffs_sorted(:, current_AC_eig_inds);
    else
        continue;
    end
end

%------------------------ Reconstruction (Geometry) ----------------------%

%Since the dimensions of eigenvectors from different octree cells may be
%different (due to different numbers of occupied voxels in different cells),
%we need to reconstruct separately for each cell.

%For each occupied cell at the chosen octree level ...
for n = 1:total_occupied_cells
    %Reconstruct the x, y, z locations of all the occupied voxels in this
    %octree cell, at the chosen level of the octree
    recon_xyz_temp = spectral_coeffs{n}*eigenvectors{n}';    %3 x N matrix   
    recon_xyz_cell{n} = recon_xyz_temp';    
end
%Convert the recon_xyz_cell cell array into an N x 3 matrix, because the
%input point cloud is in this format too
recon_xyz = cell2mat(recon_xyz_cell'); 

disp('------------------------------------------------------------');
disp('Decoder finished.');
disp('============================================================');

