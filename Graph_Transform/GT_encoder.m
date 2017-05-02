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

%p: Percentage of the largest-magnitude spectral coefficients to use for 
%   point cloud reconstruction, from all the occupied cells at octree level 
%   GT_block_lvl. NOTE: p represents the percentage of ALL spectral 
%   coefficients selected, not coefficients corresponding to x, y, and z 
%   components separately. p = 0 corresponds to no coefficients being 
%   selected for reconstruction; p = 100 corresponds to ALL coefficients 
%   being selected.

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

%practical_bits: Vector of bitrates for the transmitted data, computed by 
%                using existing entropy coders, arranged in a 2-column 
%                vector in the format: bpv total_bits. 

%entropy_bits: Vector of bitrates for the transmitted data, computed by
%              using estimates of the theoretical (Shannon) entropy, 
%              arranged in a 2-column vector in the format: bpv total_bits. 

%-------------------------------------------------------------------------%

function [nbr_occ_voxels_diff_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, practical_bits, entropy_bits] = GT_encoder(ptcloud_file, b, GT_block_lvl, p, q_stepsize, A_construction)

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

%Initialize cell arrays to hold the eigenvector and eigenvalue matrices for
%each occupied cell (node) of the octree myOT at level GT_block_lvl
eigenvectors = cell(1, total_occupied_cells);
eigenvalues = cell(1, total_occupied_cells);

%Initialize a cell array to hold the spectral coefficients for each
%occupied cell (node) of the octree myOT at level GT_block_lvl
spectral_coeffs = cell(1, total_occupied_cells);

%Read in the triangle soup vertices for the input point cloud, at a given
%octree level (level 7 in this case, where the root is considered level 1)
vertices = importdata('\\pandora\Storage\users\phil\maja\longdress1180_voxelized10_vertices.txt');

%Read in the triangle soup faces corresponding to the above vertices
faces = importdata('\\pandora\Storage\users\phil\maja\longdress1180_voxelized10_faces.txt');
%Add a 1 to each vertex index in "faces", since they start from 0 but
%MATLAB indexing starts from 1
faces = faces + 1;

%For all of the triangle faces, read in the indices of the octree cells to 
%which these faces belong
faceblocks = importdata('\\pandora\Storage\users\phil\maja\longdress1180_voxelized10_faceblocks.txt');
%Extract just the unique octree cell indices for each triangle face. The
%number of these unique indices should be equal to the number of
%total_occupied_cells.
unique_faceblocks = unique(faceblocks);

%Read in the attributes (x, y, z, r, g, b) of the refined vertices of the
%refined triangles, for the triangle soup above
attributes = importdata('\\pandora\Storage\users\phil\maja\longdress1180_voxelized10_attributes.txt');

%------------------- Individual Octree Cell Processing -------------------%

%For each occupied cell at level GT_block_lvl ...
for n = 1:total_occupied_cells
    
    %Get the total number of occupied voxels in the current cell
    nbr_occ_voxels = nbr_occ_voxels_vec(n);
        
    %------------------------ A, D, L Construction -----------------------%

    %For a fully connected graph, vertex degrees will just be equal to the
    %number of occupied voxels in the cell minus 1
    vtx_degs = (nbr_occ_voxels - 1)*ones(nbr_occ_voxels, 1); 
    
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
    
    %--------------------- Graph Transform (Geometry) --------------------%
    
    switch A_construction
        case 'fully_connected'
            %Get the set of occupied voxel coordinates (x, y, z) in the 
            %current octree cell
            current_voxel_set = occupied_voxel_coords_GT{:, n}; %N x 3 matrix
            %Project the geometry vectors (X, Y, Z) for the occupied voxels 
            %in current_voxel_set, onto the matrix of sorted eigenvectors 
            %(basis vectors) for this octree cell
            spectral_coeffs{n} = current_voxel_set'*eigvecs_sorted;  
        case 'triangle_soup'
            %Organize the vertex indices in current_faces, in ascending
            %order
            vert_inds_sorted = sort(current_faces(:));
            %Get the (x, y, z) coordinates of the set of vertices 
            %corresponding to the unique vertex indices in 
            %"vert_inds_sorted"
            current_vertices = vertices(unique(vert_inds_sorted), :);
            %Project the X, Y, and Z geometry vectors of the current vertex
            %set, onto the matrix of sorted eigenvectors (basis vectors)
            %for this octree cell
            spectral_coeffs{n} = current_vertices'*eigvecs_sorted; 
    end
end

%--------------- Spectral Coefficient Selection and Sorting --------------%
    
%Select a certain percentage (p) of the largest-magnitude spectral 
%coefficients to use for point cloud reconstruction, from ALL of the 
%coefficients that were produced for the occupied cells at octree level 
%GT_block_lvl. NOTE: p represents the percentage of ALL spectral 
%coefficients selected, not coefficients corresponding to x, y, and z 
%components separately.  
spectral_coeffs_matrix = cell2mat(spectral_coeffs);  %3 x N matrix (N is total no. of occupied voxels across all octree cells at level GT_block_lvl)
nbr_coeffs_selected = ceil((p/100)*numel(spectral_coeffs_matrix));
[~, I_coeffs] = sort(abs(spectral_coeffs_matrix(:)), 'descend');
%Initialize a matrix to store the thresholded spectral coefficients
thresh_spectral_coeffs = zeros(3, size(spectral_coeffs_matrix, 2));
%Only the selected spectral coefficients will retain their values; the
%other coefficients will be set to 0
thresh_spectral_coeffs(I_coeffs(1:nbr_coeffs_selected)) = spectral_coeffs_matrix(I_coeffs(1:nbr_coeffs_selected));   %Keep original signs (+/-) of sorted coefficients

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

%Sort the thresholded spectral coefficients in the same order as their 
%corresponding sorted eigenvalues
thresh_spectral_coeffs_sorted_temp = thresh_spectral_coeffs(:, I_eigenvalues);

%Take out all the spectral coefficients corresponding to DC components (0
%eigenvalue), for encoding separately. There should be one DC coefficient
%per occupied octree cell.
DC_spectral_coeffs_temp = thresh_spectral_coeffs_sorted_temp(:, 1:total_occupied_cells);  %3 x D vector, where D is the total no. of 0 eigenvalues
thresh_spectral_coeffs_sorted_temp(:, 1:total_occupied_cells) = [];

% %Normalize each DC coefficient by dividing its value by the square root of
% %the number of occupied voxels in the corresponding octree cell. The result 
% %is the centroid of the corresponding octree cell. 
% for i = 1:size(DC_spectral_coeffs_temp, 2)
%     DC_spectral_coeffs_temp(:, i) =  DC_spectral_coeffs_temp(:, i)./(sqrt(nbr_occ_voxels_vec(i)));   
% end
%Organize the DC spectral coefficients into a 3D x 1 vector
DC_spectral_coeffs = DC_spectral_coeffs_temp(:);
%Quantize the DC spectral coefficients using uniform scalar quantization
DC_spectral_coeffs = quantize_uniform_scalar(DC_spectral_coeffs, q_stepsize);

%Organize the remaining sorted, thresholded spectral coefficients (i.e., 
%the AC coefficients) into a 3M x 1 column vector, in the format: 
%[x1 y1 z1 ... xM yM zM]', where M is the total number of non-DC
%eigenvectors
thresh_spectral_coeffs_sorted = thresh_spectral_coeffs_sorted_temp(:);
%Quantize the AC spectral coefficients produced above, using uniform scalar
%quantization
thresh_spectral_coeffs_sorted = quantize_uniform_scalar(thresh_spectral_coeffs_sorted, q_stepsize);

%Compute the differences between successive values in nbr_occ_voxels_vec,
%and store these in a separate vector
nbr_occ_voxels_diff_vec = zeros(numel(nbr_occ_voxels_vec), 1);
for i = 1:numel(nbr_occ_voxels_vec)
    if i == 1
        nbr_occ_voxels_diff_vec(i) = nbr_occ_voxels_vec(1);
    else
        nbr_occ_voxels_diff_vec(i) = nbr_occ_voxels_vec(i) - nbr_occ_voxels_vec(i - 1);
    end
end

%---------------- Distributions of Data to be Transmitted ----------------%

%No. of occupied voxels in each occupied cell at level GT_block_lvl (these
%are all unsigned, positive integers). NOTE: This vector is not actually
%being transmitted; its distribution is just plotted here out of interest.
figure;
histogram(nbr_occ_voxels_vec);   
title({'Distribution of nbr\_occ\_voxels\_vec', '(No. of Occupied Voxels in Each Occupied Cell at the Chosen Octree Level)'});

disp('------------------------------------------------------------');
disp(' ');
disp(['MEAN (nbr_occ_voxels_vec) = ' num2str(mean(nbr_occ_voxels_vec))]);
disp(['STANDARD DEVIATION (nbr_occ_voxels_vec) = ' num2str(std(nbr_occ_voxels_vec))]);
disp(['VARIANCE (nbr_occ_voxels_vec) = ' num2str(var(nbr_occ_voxels_vec))]);
disp(' ');

%Difference vector of the values in nbr_occ_voxels_vec
figure;
histogram(nbr_occ_voxels_diff_vec);
title({'Distribution of nbr\_occ\_voxels\_diff\_vec', '(Differences between the Numbers of Occupied Voxels in', 'Successive Occupied Cells at the Chosen Octree Level)'});

disp(['MEAN (nbr_occ_voxels_diff_vec) = ' num2str(mean(nbr_occ_voxels_diff_vec))]);
disp(['STANDARD DEVIATION (nbr_occ_voxels_diff_vec) = ' num2str(std(nbr_occ_voxels_diff_vec))]);
disp(['VARIANCE (nbr_occ_voxels_diff_vec) = ' num2str(var(nbr_occ_voxels_diff_vec))]);
disp(' ');

%Thresholded, sorted (from smallest -> largest-magnitude eigenvalue) and
%quantized set of AC spectral coefficients (can be negative or non-negative
%integers)
figure;
histogram(thresh_spectral_coeffs_sorted);   
if q_stepsize == 0
    title({'Distribution of thresh\_spectral\_coeffs\_sorted', '(Thresholded, Unquantized AC Spectral Coefficients)'});
else 
    title({'Distribution of thresh\_spectral\_coeffs\_sorted', ['(Thresholded AC Spectral Coefficients, Quantized with Stepsize = ' num2str(q_stepsize) ')']});
end

disp(['MEAN (thresh_spectral_coeffs_sorted) = ' num2str(mean(thresh_spectral_coeffs_sorted))]);
disp(['STANDARD DEVIATION (thresh_spectral_coeffs_sorted) = ' num2str(std(thresh_spectral_coeffs_sorted))]);
disp(['VARIANCE (thresh_spectral_coeffs_sorted) = ' num2str(var(thresh_spectral_coeffs_sorted))]);
disp(' ');

%Quantized DC spectral coefficients
figure;
histogram(DC_spectral_coeffs);
if q_stepsize == 0
    title({'Distribution of DC\_spectral\_coeffs', '(Unquantized DC Spectral Coefficients)'});
else
    title({'Distribution of DC\_spectral\_coeffs', ['(DC Spectral Coefficients, Quantized with Stepsize = ' num2str(q_stepsize) ')']});
end

disp(['MEAN (DC_spectral_coeffs) = ' num2str(mean(DC_spectral_coeffs))]);
disp(['STANDARD DEVIATION (DC_spectral_coeffs) = ' num2str(std(DC_spectral_coeffs))]);
disp(['VARIANCE (DC_spectral_coeffs) = ' num2str(var(DC_spectral_coeffs))]);
disp(' ');

%-------------------- Encoding Data to be Transmitted --------------------%

%NOTES: rlgr() is Adaptive Run-Length Golomb Rice encoding; entropy()
%computes an estimate of the theoretical (practical) Shannon entropy of a
%given sequence of symbols. 

%Encode the data to be transmitted to the decoder (i.e., figure out how
%many bits would be needed to encode this data) ...

%Differences between the numbers of occupied voxels in successive occupied
%cells at level GT_block_lvl (can be negative or non-negative integers)
bits_nbr_occ_voxels_diff_rlgr = rlgr(nbr_occ_voxels_diff_vec);    
disp('------------------------------------------------------------');
disp(['bits_nbr_occ_voxels_diff_rlgr: ' num2str(bits_nbr_occ_voxels_diff_rlgr)]);
disp(['bpv: ' num2str(bits_nbr_occ_voxels_diff_rlgr/sum(nbr_occ_voxels_vec))]);

min_practical_bits_nbr_occ_voxels_diff = min(bits_nbr_occ_voxels_diff_rlgr);
disp(['min_practical_bits_nbr_occ_voxels_diff: ' num2str(min_practical_bits_nbr_occ_voxels_diff)]);
disp(['bpv: ' num2str(min_practical_bits_nbr_occ_voxels_diff/sum(nbr_occ_voxels_vec))]);

bits_nbr_occ_voxels_diff_entropy_persymbol = entropy_calc(nbr_occ_voxels_diff_vec);   %Avg. no. of bits per symbol
bits_nbr_occ_voxels_diff_entropy = bits_nbr_occ_voxels_diff_entropy_persymbol*length(nbr_occ_voxels_diff_vec); %Total no. of bits for all symbols
disp(['bits_nbr_occ_voxels_diff_entropy: ' num2str(bits_nbr_occ_voxels_diff_entropy) ' (' num2str(bits_nbr_occ_voxels_diff_entropy_persymbol) ' bits per symbol)']);
disp(['bpv: ' num2str(bits_nbr_occ_voxels_diff_entropy/sum(nbr_occ_voxels_vec))]);
disp(' ');

%Quantized, thresholded and sorted (from smallest -> largest-magnitude 
%eigenvalue) set of AC spectral coefficients (can be negative or 
%non-negative integers)
bits_quant_coeffs_rlgr = rlgr(thresh_spectral_coeffs_sorted);
disp(['bits_quant_coeffs_rlgr: ' num2str(bits_quant_coeffs_rlgr)]);
disp(['bpv: ' num2str(bits_quant_coeffs_rlgr/sum(nbr_occ_voxels_vec))]);

bits_quant_coeffs_entropy_persymbol = entropy_calc(thresh_spectral_coeffs_sorted);  %Avg. no. of bits per symbol
bits_quant_coeffs_entropy = bits_quant_coeffs_entropy_persymbol*length(thresh_spectral_coeffs_sorted);  %Total no. of bits for all symbols
disp(['bits_quant_coeffs_entropy: ' num2str(bits_quant_coeffs_entropy) ' (' num2str(bits_quant_coeffs_entropy_persymbol) ' bits per symbol)']);
disp(['bpv: ' num2str(bits_quant_coeffs_entropy/sum(nbr_occ_voxels_vec))]);
disp(' ');

%Quantized DC spectral coefficients
bits_DC_coeffs_rlgr = rlgr(DC_spectral_coeffs);
disp(['bits_DC_coeffs_rlgr: ' num2str(bits_DC_coeffs_rlgr)]);
disp(['bpv: ' num2str(bits_DC_coeffs_rlgr/sum(nbr_occ_voxels_vec))]);

bits_DC_coeffs_entropy_persymbol = entropy_calc(DC_spectral_coeffs);    %Avg. no. of bits per symbol
bits_DC_coeffs_entropy = bits_DC_coeffs_entropy_persymbol*length(DC_spectral_coeffs);   %Total no. of bits for all symbols
disp(['bits_DC_coeffs_entropy: ' num2str(bits_DC_coeffs_entropy) ' (' num2str(bits_DC_coeffs_entropy_persymbol) ' bits per symbol)']);
disp(['bpv: ' num2str(bits_DC_coeffs_entropy/sum(nbr_occ_voxels_vec))]);
disp(' ');

%TOTAL BITS
total_bits_practical = min_practical_bits_nbr_occ_voxels_diff + bits_quant_coeffs_rlgr + bits_DC_coeffs_rlgr; 
disp('------------------------------------------------------------');
disp(['TOTAL BITS (PRACTICAL): ' num2str(total_bits_practical)]);
total_bits_entropy = bits_nbr_occ_voxels_diff_entropy + bits_quant_coeffs_entropy + bits_DC_coeffs_entropy; 
disp(['TOTAL BITS (ENTROPY): ' num2str(total_bits_entropy)]);
disp(' ');

%Bits per voxel
bpv_practical = total_bits_practical/sum(nbr_occ_voxels_vec);
disp(['bpv_practical: ' num2str(bpv_practical) ' (' num2str(100*((min_practical_bits_nbr_occ_voxels_diff/sum(nbr_occ_voxels_vec))/bpv_practical)) '% for nbr. occ. vox. diffs, ' num2str(100*((bits_quant_coeffs_rlgr/sum(nbr_occ_voxels_vec))/bpv_practical)) '% for AC coeffs, ' num2str(100*((bits_DC_coeffs_rlgr/sum(nbr_occ_voxels_vec))/bpv_practical)) '% for DC coeffs)']);
bpv_entropy = total_bits_entropy/sum(nbr_occ_voxels_vec);
disp(['bpv_entropy: ' num2str(bpv_entropy) ' (' num2str(100*((bits_nbr_occ_voxels_diff_entropy/sum(nbr_occ_voxels_vec))/bpv_entropy)) '% for nbr. occ. vox. diffs, ' num2str(100*((bits_quant_coeffs_entropy/sum(nbr_occ_voxels_vec))/bpv_entropy)) '% for AC coeffs, ' num2str(100*((bits_DC_coeffs_entropy/sum(nbr_occ_voxels_vec))/bpv_entropy)) '% for DC coeffs)']);

%Store the bpv values and the total_bits values in their respective vectors
practical_bits = [bpv_practical, total_bits_practical];
entropy_bits = [bpv_entropy, total_bits_entropy];

disp('------------------------------------------------------------');
disp('Encoder finished.');
disp('============================================================');




