%Script for running the GT codec with different parameters, to obtain 
%reconstructions of a given input point cloud at different bitrates.

%----------------------------- User Inputs -------------------------------%

%Path to where the input point cloud files are contained (use backslashes, 
%including a backslash at the end of the path).
data_path = '\\Pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\';
%Path to the directory that contains the codec outputs (use backslashes, 
%including a backslash at the end of the path).
codec_results_path = '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\';
%Name of the input point cloud that you wish to use (don't include the 
%_voxelizedN or .ply file extension in the name).
ptcloud_name = 'boxer'; %Must be in PLY format for now
%The voxelization level for the input point cloud. If no voxelization, put 
%0 here. Write the voxelizedN number as a string, e.g., '10'.
voxelizedN = '10';
%Bit depth for Morton codes and octree. b also determines the number of 
%levels in the octree that will be generated (apart from the root level). 
%The total number of octree levels INCLUDING the root level will therefore 
%be: b + 1. IMPORTANT: If using a voxelized point cloud as input, b must be 
%equal to the voxelization level used to produce that point cloud, e.g., 
%b = 10 for voxelized10 point clouds, b = 11 for voxelized11 clouds, etc.
b = 10;
%Desired octree level, at which the Graph Transform will be computed for 
%each occupied cell. IMPORTANT: GT_block_lvl must be <= b.
GT_block_lvl = 8;
%Flag to indicate whether to use a fixed p value (see below) or a fixed
%q_stepsize (see below). fixed_p_qstep = 0 corresponds to using a fixed
%p value (and varying q_stepsize values); fixed_p_qstep = 1 corresponds to
%using a fixed q_stepsize (and varying p values).
fixed_p_qstep = 1;
%Percentage of the largest spectral coefficients to use for point cloud 
%reconstruction, from each occupied cell at octree level GT_block_lvl, for 
%which the Graph Transform is computed. NOTE: p represents the percentage 
%of ALL spectral coefficients selected, not coefficients corresponding to 
%x, y, and z separately. p = 0 corresponds to no coefficients being 
%selected for reconstruction; p = 100 corresponds to ALL coefficients being
%selected.
p = 100;    %Comment this out if using varying p values
%Step size used for uniform scalar quantization. Values for step_size 
%should ideally be powers of 2. q_stepsize = 1 corresponds to just rounding 
%the input to the nearest integer, so the least amount of (uniform)
%quantization.  
q_stepsize = 2; %Comment this out if using various quantization step sizes

%-------------------------------------------------------------------------%

if str2double(voxelizedN) == 0
    vox_novox_dir = 'OriginalPLY';
    reconstruction_prefix = ptcloud_name;
else
    vox_novox_dir = ['voxelized' voxelizedN];
    reconstruction_prefix = [ptcloud_name '_voxelized' voxelizedN];
end

ptcloud_file = [data_path vox_novox_dir '\' reconstruction_prefix '.ply'];

%Read in input point cloud, so that we can extract properties of the PLY
%file, which will be needed later
[plyStruct, A, format] = plyRead(ptcloud_file);
%Make a copy of plyStruct, as it will be modified
plyStruct2 = plyStruct;

%Initialize counter that will keep track of how many point cloud
%reconstructions we have done so far
reconstruction_counter = 1;

%If using a fixed q_stepsize
if (fixed_p_qstep == 1)
    %Define the codec names, which will be used to write the results to the
    %correct directories
    codec_name_practical = ['GTGeom_QS' num2str(q_stepsize) '_B' num2str(GT_block_lvl) '_PRACTICAL'];
    codec_name_entropy = ['GTGeom_QS' num2str(q_stepsize) '_B' num2str(GT_block_lvl) '_ENTROPY'];
    %Define the output directories for each codec
    outputdir_practical = [codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name_practical '\'];
    outputdir_entropy = [codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name_entropy '\'];
    %Create the directories corresponding to the above codec names
    %(NOTE: if these directories already exist, the operation below
    %will not overwrite the contents inside them; MATLAB will just
    %generate a warning saying that these directories already exist)
    newdir_stat_practical = mkdir(outputdir_practical); 
    if newdir_stat_practical == 1
        disp(['Created output directory: ' outputdir_practical]);
    else
        error(['ERROR: Could not create output directory ' outputdir_practical]);
    end
    newdir_stat_entropy = mkdir(outputdir_entropy); 
    if newdir_stat_entropy == 1
        disp(['Created output directory: ' outputdir_entropy]);
    else
        error(['ERROR: Could not create output directory ' outputdir_entropy]);
    end
    
    %Open a text file in each of the directories created above, to write 
    %the bitrates into. If these text files already exist, delete them and 
    %open new ones.
    
    %For practical bitrates
    fid_practical = fopen([outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);
    if fid_practical ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_practical);
        %Delete the existing file
        disp(['The file ' outputdir_practical reconstruction_prefix '_geom_bitrates.txt already exists. Deleting ...']);
        delete([outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_practical = fopen([outputdir_practical reconstruction_prefix '_geom_bitrates.txt'], 'a');
    if fid_practical == -1
        error(['ERROR: Could not open text file ' outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);
    else
        disp(['Successfully opened file ' outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);
    end

    %For entropy bitrates
    fid_entropy = fopen([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);
    if fid_entropy ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_entropy);
        %Delete the existing file
        disp(['The file ' outputdir_entropy reconstruction_prefix '_geom_bitrates.txt already exists. Deleting ...']);
        delete([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_entropy = fopen([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt'], 'a');
    if fid_entropy == -1
        error(['ERROR: Could not open text file ' outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);
    else
        disp(['Successfully opened file ' outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);
    end
    
    %Vary the p values used for spectral coefficient selection, to obtain
    %the different rates on the rate-distortion plots (NOTE: Remember that
    %the highest bitrate, therefore the best reconstruction, should come 
    %first)
    for p = [50 40 30 20 10 5]
        disp('------------------------------------------------------------');
        disp(['Processing p = ' num2str(p) ' ...']);
        disp('------------------------------------------------------------');
        %Run the encoder
        [nbr_occ_voxels_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, practical_bits, entropy_bits] = GT_encoder(ptcloud_file, b, GT_block_lvl, p, q_stepsize);        
        %Write the bpv (bits per reconstructed voxel, or point) and total     
        %bitrate to file, for the current compressed point cloud
        fprintf(fid_practical, '%f %f\r\n', practical_bits);
        fprintf(fid_entropy, '%f %f\r\n', entropy_bits);
        %Run the decoder
        recon_xyz = GT_decoder(nbr_occ_voxels_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, q_stepsize);
        plyStruct2.propArrayListList{1}{1} = recon_xyz(:, 1);   %X coordinates
        plyStruct2.propArrayListList{1}{2} = recon_xyz(:, 2);   %Y coordinates
        plyStruct2.propArrayListList{1}{3} = recon_xyz(:, 3);   %Z coordinates
        %Write the current reconstructed point cloud to file, in PLY format
        if (reconstruction_counter <= 9)
            plyWrite(plyStruct2, [outputdir_practical reconstruction_prefix '_distorted0' num2str(reconstruction_counter) '.ply'], format);
            plyWrite(plyStruct2, [outputdir_entropy reconstruction_prefix '_distorted0' num2str(reconstruction_counter) '.ply'], format);
        else
            plyWrite(plyStruct2, [outputdir_practical reconstruction_prefix '_distorted' num2str(reconstruction_counter) '.ply'], format);
            plyWrite(plyStruct2, [outputdir_entropy reconstruction_prefix '_distorted' num2str(reconstruction_counter) '.ply'], format);
        end
        reconstruction_counter = reconstruction_counter + 1;
    end
    %Close the open bitrates text files
    fclose(fid_practical);
    fclose(fid_entropy);
    
%If using a fixed p value
elseif (fixed_p_qstep == 0)
    %Define the codec names, which will be used to write the results to
    %the current directories
    codec_name_practical = ['GTGeom_p' num2str(p) '_B' num2str(GT_block_lvl) '_PRACTICAL'];
    codec_name_entropy = ['GTGeom_p' num2str(p) '_B' num2str(GT_block_lvl) '_ENTROPY']; 
    %Define the output directories for each codec
    outputdir_practical = [codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name_practical '\'];
    outputdir_entropy = [codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name_entropy '\'];
    %Create the directories corresponding to the above codec names
    %(NOTE: if these directories already exist, the operation below
    %will not overwrite the contents inside them; MATLAB will just
    %generate a warning saying that these directories already exist)
    newdir_stat_practical = mkdir(outputdir_practical); 
    if newdir_stat_practical == 1
        disp(['Created output directory: ' outputdir_practical]);
    else
        disp(['ERROR: Could not create output directory ' outputdir_practical]);
    end
    newdir_stat_entropy = mkdir(outputdir_entropy); 
    if newdir_stat_entropy == 1
        disp(['Created output directory: ' outputdir_entropy]);
    else
        disp(['ERROR: Could not create output directory ' outputdir_entropy]);
    end
    
    %Open a text file in each of the directories created above, to write 
    %the bitrates into. If these text files already exist, delete them and 
    %open new ones.
    
    %For practical bitrates
    fid_practical = fopen([outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);
    if fid_practical ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_practical);
        %Delete the existing file
        disp(['The file ' outputdir_practical reconstruction_prefix '_geom_bitrates.txt already exists. Deleting ...']);
        delete([outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_practical = fopen([outputdir_practical reconstruction_prefix '_geom_bitrates.txt'], 'a');
    if fid_practical == -1
        error(['ERROR: Could not open text file ' outputdir_practical reconstruction_prefix '_geom_bitrates.txt']);
    end

    %For entropy bitrates
    fid_entropy = fopen([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);
    if fid_entropy ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_entropy);
        %Delete the existing file
        disp(['The file ' outputdir_entropy reconstruction_prefix '_geom_bitrates.txt already exists. Deleting ...']);
        delete([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_entropy = fopen([outputdir_entropy reconstruction_prefix '_geom_bitrates.txt'], 'a');
    if fid_entropy == -1
        error(['ERROR: Could not open text file ' outputdir_entropy reconstruction_prefix '_geom_bitrates.txt']);
    end
   
    %Vary the q_stepsize values used for spectral coefficient selection, to 
    %obtain the different rates on the rate-distortion plots (NOTE:
    %Remember that the highest bitrate, therefore the best reconstruction,
    %should come first)
    for q_stepsize = [1 2 4 8 16 32 64 128]
        disp('------------------------------------------------------------');
        disp(['Processing q_stepsize = ' num2str(q_stepsize) ' ...']);
        disp('------------------------------------------------------------');
        %Run the encoder
        [nbr_occ_voxels_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, practical_bits, entropy_bits] = GT_encoder(ptcloud_file, b, GT_block_lvl, p, q_stepsize);        
        %Write the bpv (bits per reconstructed voxel, or point) and total 
        %bitrate to file, for the current compressed point cloud
        fprintf(fid_practical, '%f %f\r\n', practical_bits);
        fprintf(fid_entropy, '%f %f\r\n', entropy_bits);
        %Run the decoder
        recon_xyz = GT_decoder(nbr_occ_voxels_vec, thresh_spectral_coeffs_sorted, DC_spectral_coeffs, q_stepsize);
        plyStruct2.propArrayListList{1}{1} = recon_xyz(:, 1);   %X coordinates
        plyStruct2.propArrayListList{1}{2} = recon_xyz(:, 2);   %Y coordinates
        plyStruct2.propArrayListList{1}{3} = recon_xyz(:, 3);   %Z coordinates
        %Write the current reconstructed point cloud to file, in PLY format
        if (reconstruction_counter <= 9)
            plyWrite(plyStruct2, [outputdir_practical reconstruction_prefix '_distorted0' num2str(reconstruction_counter) '.ply'], format);
            plyWrite(plyStruct2, [outputdir_entropy reconstruction_prefix '_distorted0' num2str(reconstruction_counter) '.ply'], format);
        else
            plyWrite(plyStruct2, [outputdir_practical reconstruction_prefix '_distorted' num2str(reconstruction_counter) '.ply'], format);
            plyWrite(plyStruct2, [outputdir_entropy reconstruction_prefix '_distorted' num2str(reconstruction_counter) '.ply'], format);
        end
        reconstruction_counter = reconstruction_counter + 1;        
    end
    %Close the open bitrates text files
    fclose(fid_practical);
    fclose(fid_entropy);
end


    

    
    
    
    