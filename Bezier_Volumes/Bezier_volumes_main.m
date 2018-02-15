%Need the following directories (add these directories and all of their
%sub-directories to the current MATLAB path): Bezier_Volumes (run this
%script file from inside this directory), Encoding, Octree, Phil, 
%Quantization.

%------------------------------ User Inputs ------------------------------%

%Name of the input point cloud (do not include the _voxelizedN or .ply file 
%extension in the name)
ptcloud_name = 'longdress_1300'; %Must be in PLY format
%Bit depth for Morton codes and octree. b also determines the number of 
%levels in the octree that will be generated (apart from the root level). 
%The total number of octree levels INCLUDING the root level will therefore 
%be b + 1. IMPORTANT: If using a voxelized point cloud as input, b must be 
%equal to the voxelization level used to produce that point cloud, 
%e.g., b = 10 for voxelized10 point clouds, b = 11 for voxelized11 clouds,
%etc.
b = 10;
%Full path to the input PLY file
%ptcloud_file = ['\\pandora\storage\users\phil\maja\voxelized7_Test\' ptcloud_name '_voxelized' num2str(b) '.ply'];
ptcloud_file = ['\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized' num2str(b) '_WithNormalsAndCentroids\' ptcloud_name '_voxelized' num2str(b) '.ply'];
%ptcloud_file = ['\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized' num2str(b) '_WithNormals\' ptcloud_name '_voxelized' num2str(b) '.ply'];
%Octree level(s) to use (one at a time, if there is more than one listed
%below) as the base level for transmitting control points, and from which
%the wavelet analysis will start. start_OT_lvl can go from the root level
%(start_OT_lvl = 1) up to 2 levels before the leaf level
%(start_OT_lvl = b - 1)), or up to (max_lvl - 1) if max_lvl is not the leaf 
%level. But because the lowest octree levels usually do not contain zero 
%crossings (since the octree cells here are quite large and so all the cell 
%corners are outside the point cloud), start_OT_lvl should be adjusted to 
%start at the lowest level that contains zero crossings. 
start_OT_lvl = 3;
%Highest octree level at which the Bezier control points should be computed
%at the encoder and for which the wavelet coefficients should be sent to 
%the decoder. max_OT_lvl must go from (start_OT_lvl + 1) and can go up to 
%b + 1. 
max_OT_lvl = b + 1;  %Write numbers in DEscending order, because file ..._distorted01.ply must correspond to the best reconstruction (and highest bitrate)
%Quantization stepsize for uniform scalar quantization of the control 
%points at the chosen base level (start_lvl) and for all of the wavelet
%coefficients that will be computed at the encoder
%q_stepsizes = [0.25 0.5 1 2 3];
q_stepsizes = 5;
%Decide whether or not to prune the octree cells at the encoder, which
%contain zero wavelet coefficients on all of their corners, and therefore
%whether to prune the corresponding wavelet coefficient tree: 
%prune_flag = 1 => prune; prune_flag = 0 => do not prune
prune_flag = 1;
prune_levels_all = [7 8 9];
%Define threshold for pruning quantized wavelet coefficients
%(zero_threshold can only have integer values, because the quantized
%wavelet coefficients will all be integers)
%zero_threshold_for_pruning = [0 1 2 3 4];
zero_threshold_for_pruning = 0;
%Choose whether to run in "debug" mode (debug_flag = 1), where a lot of 
%information is written to the log file and certain figures are displayed,
%which can be used to check that the program is working correctly. If
%debug_flag = 0, only a few basic outputs are written to log file, which 
%also speeds up execution of the code.
debug_flag = 1;
%Choose whether to apply colour compression on top of the BV-compressed
%geometry (colour_compression = 1), or not (colour_compression = 0)
colour_compression = 0;
%Name of the colour codec that will be used (if colour_compression = 1)
col_codec = 'RAHT';

%-------------------------------------------------------------------------%

%Create the output directory where the reconstructed point cloud(s) from
%this codec will be written
if colour_compression == 0
    output_dir = ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\'];
else
    output_dir = ['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume_' col_codec '\'];
end
newdir_stat = mkdir(output_dir);
if newdir_stat == 1
    disp(' ');
    disp(['Created PLY output directory: ' output_dir]);
else
    disp(' ');
    error(['ERROR: Could not create output directory ' output_dir]);
end
disp(' ');

%Start recording a log of all MATLAB command line input and output from
%this point onwards. But first check if a log file already exists and
%delete it if it does, then create a new log file.
fid_log = fopen(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\log_BezVol.txt']);
if fid_log ~= -1
    %Close the open file (must close before deleting)
    fclose(fid_log);
    %Delete the existing file
    disp(['The file ' '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\log_BezVol.txt already exists. Deleting ...']);
    delete(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\log_BezVol.txt']); 
end
%Create a new log file
diary(['\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\log_BezVol.txt']);
disp(['Created and opened text file for writing: \\Pandora\builds\test\Data\Compression\PLY\Codec_Results\' ptcloud_name '\voxelized' num2str(b) '\BezierVolume\log_BezVol.txt']);
disp('------------------------------------------------------------');
%Start recording command line input and output
diary on;
%Write in the user inputs first
disp(['debug_flag = ' num2str(debug_flag)]);
disp(' ');
disp(['ptcloud_name = ' ptcloud_name]);
disp(['b = ' num2str(b)]);
disp(['ptcloud_file = ' ptcloud_file]);
%disp(['max_lvl = ' num2str(max_lvl)]);
%disp(['q_stepsize = ' num2str(q_stepsize)]);
disp('------------------------------------------------------------');

%Open a text file inside output_dir, to write the geometry bitrates into.
%If this text file already exists, delete it and create a new one.
fid_geom_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
if fid_geom_bits ~= -1
    %Close the open file (must close before deleting)
    fclose(fid_geom_bits);
    %Delete the existing file
    disp(['The file ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt already exists. Deleting ...']);
    delete([output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);  
end
%Open a new text file in append mode
fid_geom_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt'], 'a');
if fid_geom_bits == -1
    error(['ERROR: Could not open text file ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
else
    disp(['Created and opened text file for writing: ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
end
disp('------------------------------------------------------------');

if colour_compression == 1
    %Open a text file inside output_dir, to write the colour bitrates into.
    %If this text file already exists, delete it and create a new one.
    fid_col_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);
    if fid_col_bits ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_col_bits);
        %Delete the existing file
        disp(['The file ' output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt already exists. Deleting ...']);
        delete([output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_col_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt'], 'a');
    if fid_col_bits == -1
        error(['ERROR: Could not open text file ' output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);
    else
        disp(['Created and opened text file for writing: ' output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);
    end
    disp('------------------------------------------------------------');

    %Open a text file inside output_dir, to write the total (geom. + col.)
    %bitrates into. If this text file already exists, delete it and create 
    %a new one.
    fid_all_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);
    if fid_all_bits ~= -1
        %Close the open file (must close before deleting)
        fclose(fid_all_bits);
        %Delete the existing file
        disp(['The file ' output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt already exists. Deleting ...']);
        delete([output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);  
    end
    %Open a new text file in append mode
    fid_all_bits = fopen([output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt'], 'a');
    if fid_all_bits == -1
        error(['ERROR: Could not open text file ' output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);
    else
        disp(['Created and opened text file for writing: ' output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);
    end
    disp('------------------------------------------------------------');
end %End check if colour_compression == 1

%Read in the input point cloud, so that we can extract properties of
%the PLY file
[plyStruct, A, format] = plyRead(ptcloud_file);
%Make a copy of plyStruct, as it will be modified
plyStruct2 = plyStruct;

%If using variable-level octree pruning, use the different values specified
%in zero_threshold, to obtain different reconstruction points for a 
%rate-distortion curve
if (prune_flag == 1) && (~isempty(zero_threshold_for_pruning))
    start_lvl = start_OT_lvl;
    max_lvl = max_OT_lvl;
    %zero_threshold_cntr = 1;
    q_stepsize_cntr = 1;
    %To get separate R-D curves for each zero_threshold value (but using
    %all q_stepsize values to get the operation points in each curve), use
    %one zero_threshold value at a time below, and just vary q_stepsize
    for zero_threshold = zero_threshold_for_pruning
        disp(' ');
        disp(['prune_flag = ' num2str(prune_flag)]);
        disp(['zero_threshold = ' num2str(zero_threshold)]);
        for q_stepsize = q_stepsizes
            disp(['q_stepsize = ' num2str(q_stepsize)]);
            
            %Run encoder
            [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points, post_pruning_array_forDec] = Bezier_volumes_encoder(debug_flag, ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag, colour_compression, zero_threshold);
            
            %Run decoder
            [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(debug_flag, occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, prune_flag, colour_compression, post_pruning_array_forDec);

            %Create a new cell array for the plyStruct2 property arrays, 
            %which contains the reconstructed voxel x, y, z coordinates
            %and, if colour compression was used, their corresponding 
            %decoded colours (R, G, B) 
            plyStruct2.propArrayListList = cell(1, 1);
            plyStruct2.propArrayListList{1}{1} = reconstructed_vox_pos(:, 1);   %Reconstructed voxel X coordinates
            plyStruct2.propArrayListList{1}{2} = reconstructed_vox_pos(:, 2);   %Reconstructed voxel Y coordinates
            plyStruct2.propArrayListList{1}{3} = reconstructed_vox_pos(:, 3);   %Reconstructed voxel Z coordinates
            if colour_compression == 1
                plyStruct2.propArrayListList{1}{4} = reconstructed_vox_pos(:, 4);   %R colour value
                plyStruct2.propArrayListList{1}{5} = reconstructed_vox_pos(:, 5);   %G colour value
                plyStruct2.propArrayListList{1}{6} = reconstructed_vox_pos(:, 6);   %B colour value
            end
            %Create a new cell array for the plyStruct2 property types, 
            %which contains the data types of the reconstructed voxel x, y,
            %z coordinates and, if colour compression was used, the data
            %types of the colour data 
            plyStruct2.propTypeListList = cell(1, 1);
            plyStruct2.propTypeListList{1}(1) = "float";
            plyStruct2.propTypeListList{1}(2) = "float";
            plyStruct2.propTypeListList{1}(3) = "float";
            if colour_compression == 1
                plyStruct2.propTypeListList{1}(4) = "uchar";
                plyStruct2.propTypeListList{1}(5) = "uchar";
                plyStruct2.propTypeListList{1}(6) = "uchar";        
            end
            %Create a new cell array for the plyStruct2 property names, 
            %which contains the names for the reconstructed voxel x, y, z 
            %coordinates and, if colour compression was used, the colour
            %data 
            plyStruct2.propNameListList = cell(1, 1);
            plyStruct2.propNameListList{1}(1) = "x";
            plyStruct2.propNameListList{1}(2) = "y";
            plyStruct2.propNameListList{1}(3) = "z";
            if colour_compression == 1
                plyStruct2.propNameListList{1}(4) = "red";
                plyStruct2.propNameListList{1}(5) = "green";
                plyStruct2.propNameListList{1}(6) = "blue";         
            end

            if q_stepsize_cntr <= 9
                plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(q_stepsize_cntr) '.ply'], format);
                disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(q_stepsize_cntr) '.ply']);
            else
                plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(q_stepsize_cntr) '.ply'], format);
                disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(q_stepsize_cntr) '.ply']);
            end
            disp(' ');

            %Write the geometry bitrates to text file (the highest bitrate 
            %should come first, corresponding to the best reconstruction in
            %..._distorted01.ply)
            fprintf(fid_geom_bits, '%f %f\r\n', [total_geom_bpv total_geom_bits]);
            disp(['Finished writing to geometry bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
            disp('------------------------------------------------------------');

            if colour_compression == 1
                %Write the colour bitrates to text file (the highest 
                %bitrate should come first, corresponding to the best 
                %reconstruction in ..._distorted01.ply)
                fprintf(fid_col_bits, '%f %f\r\n', [total_col_bpv total_col_bits]);
                disp(['Finished writing to colour bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);
                disp('------------------------------------------------------------');            

                %Write the total (geometry + colour) bitrates to text file 
                %(the highest bitrate should come first, corresponding to 
                %the best reconstruction in ..._distorted01.ply)
                fprintf(fid_all_bits, '%f %f\r\n', [total_bpv total_bits]);
                disp(['Finished writing to total bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);
                disp('------------------------------------------------------------');  
            end
            
            %Close all open figures
            close all;
            
            %Update counter for the next q_stepsize
            q_stepsize_cntr =  q_stepsize_cntr + 1;
            
        end %End q_stepsize

        %Update counter for the next zero_threshold
        %zero_threshold_cntr = zero_threshold_cntr + 1;
        
    end %End zero_threshold
end %End check if (prune_flag == 1) && (~isempty(zero_threshold))

%If using a fixed-level octree pruning, use one fixed level at a time,
%while varying q_stepsize, to get different R-D curves (from different
%levels), where the different R-D points on each curve represent different
%q_stepsize values
if (prune_flag == 0)
    start_lvl = start_OT_lvl;
    max_lvl = max_OT_lvl;
    q_stepsize_cntr = 1; 
    for prune_level = prune_levels_all
        disp(' ');
        disp(['prune_flag = ' num2str(prune_flag)]);
        disp(['prune_level = ' num2str(prune_level)]);
        for q_stepsize = q_stepsizes
            disp(['q_stepsize = ' num2str(q_stepsize)]);
            
            %Run encoder
            [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points, ~] = Bezier_volumes_encoder(debug_flag, ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag, colour_compression, prune_level);
            
            %Run decoder
            [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(debug_flag, occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, prune_flag, colour_compression, prune_level);

            %Create a new cell array for the plyStruct2 property arrays, 
            %which contains the reconstructed voxel x, y, z coordinates
            %and, if colour compression was used, their corresponding 
            %decoded colours (R, G, B) 
            plyStruct2.propArrayListList = cell(1, 1);
            plyStruct2.propArrayListList{1}{1} = reconstructed_vox_pos(:, 1);   %Reconstructed voxel X coordinates
            plyStruct2.propArrayListList{1}{2} = reconstructed_vox_pos(:, 2);   %Reconstructed voxel Y coordinates
            plyStruct2.propArrayListList{1}{3} = reconstructed_vox_pos(:, 3);   %Reconstructed voxel Z coordinates
            if colour_compression == 1
                plyStruct2.propArrayListList{1}{4} = reconstructed_vox_pos(:, 4);   %R colour value
                plyStruct2.propArrayListList{1}{5} = reconstructed_vox_pos(:, 5);   %G colour value
                plyStruct2.propArrayListList{1}{6} = reconstructed_vox_pos(:, 6);   %B colour value
            end
            %Create a new cell array for the plyStruct2 property types, 
            %which contains the data types of the reconstructed voxel x, y,
            %z coordinates and, if colour compression was used, the data
            %types of the colour data 
            plyStruct2.propTypeListList = cell(1, 1);
            plyStruct2.propTypeListList{1}(1) = "float";
            plyStruct2.propTypeListList{1}(2) = "float";
            plyStruct2.propTypeListList{1}(3) = "float";
            if colour_compression == 1
                plyStruct2.propTypeListList{1}(4) = "uchar";
                plyStruct2.propTypeListList{1}(5) = "uchar";
                plyStruct2.propTypeListList{1}(6) = "uchar";        
            end
            %Create a new cell array for the plyStruct2 property names, 
            %which contains the names for the reconstructed voxel x, y, z 
            %coordinates and, if colour compression was used, the colour
            %data 
            plyStruct2.propNameListList = cell(1, 1);
            plyStruct2.propNameListList{1}(1) = "x";
            plyStruct2.propNameListList{1}(2) = "y";
            plyStruct2.propNameListList{1}(3) = "z";
            if colour_compression == 1
                plyStruct2.propNameListList{1}(4) = "red";
                plyStruct2.propNameListList{1}(5) = "green";
                plyStruct2.propNameListList{1}(6) = "blue";         
            end

            if q_stepsize_cntr <= 9
                plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(q_stepsize_cntr) '.ply'], format);
                disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(q_stepsize_cntr) '.ply']);
            else
                plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(q_stepsize_cntr) '.ply'], format);
                disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(q_stepsize_cntr) '.ply']);
            end
            disp(' ');

            %Write the geometry bitrates to text file (the highest bitrate 
            %should come first, corresponding to the best reconstruction in
            %..._distorted01.ply)
            fprintf(fid_geom_bits, '%f %f\r\n', [total_geom_bpv total_geom_bits]);
            disp(['Finished writing to geometry bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
            disp('------------------------------------------------------------');

            if colour_compression == 1
                %Write the colour bitrates to text file (the highest 
                %bitrate should come first, corresponding to the best 
                %reconstruction in ..._distorted01.ply)
                fprintf(fid_col_bits, '%f %f\r\n', [total_col_bpv total_col_bits]);
                disp(['Finished writing to colour bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_col_bitrates.txt']);
                disp('------------------------------------------------------------');            

                %Write the total (geometry + colour) bitrates to text file 
                %(the highest bitrate should come first, corresponding to 
                %the best reconstruction in ..._distorted01.ply)
                fprintf(fid_all_bits, '%f %f\r\n', [total_bpv total_bits]);
                disp(['Finished writing to total bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_all_bitrates.txt']);
                disp('------------------------------------------------------------');  
            end
            
            %Close all open figures
            close all;
            
            %Update counter for the next q_stepsize
            q_stepsize_cntr =  q_stepsize_cntr + 1;
            
        end %End q_stepsize          
    end %End prune_level
end %End check if (prune_flag == 1)

% %Try using different start levels, to obtain different reconstruction
% %points for a rate-distortion curve
% %start_lvl_cntr = 1;
% for start_lvl = start_OT_lvl
%     disp(' ');
%     disp(['start_lvl = ' num2str(start_lvl)]);
%     %Try using different end levels for the wavelet coefficients that are
%     %computed and transmitted to the decoder, to get different 
%     %reconstruction points for a rate-distortion curve
%     %max_lvl_cntr = 1;
%     for max_lvl = max_OT_lvl
%         disp(' ');
%         disp(['max_lvl = ' num2str(max_lvl)]);
%     
%         %Run encoder
%         %[occupancy_codes_forDec, post_pruning_array_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points] = Bezier_volumes_encoder(ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name);
%         if prune_flag == 1
%             [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points, post_pruning_array_forDec] = Bezier_volumes_encoder(ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag, zero_threshold);
%         else
%             [occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, total_geom_bits, total_geom_bpv, reconstructed_control_points] = Bezier_volumes_encoder(ptcloud_file, b, start_lvl, max_lvl, q_stepsize, ptcloud_name, prune_flag);
%         end
%         
%         %Run decoder
%         %[reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(occupancy_codes_forDec, post_pruning_array_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, reconstructed_control_points);
%         if prune_flag == 1
%             [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, reconstructed_control_points, prune_flag, post_pruning_array_forDec);
%         else
%             [reconstruction_decoder, reconstructed_vox_pos] = Bezier_volumes_decoder(occupancy_codes_forDec, rec_ctrlpts_forDec, wavelet_coeffs_forDec, start_lvl, max_lvl, q_stepsize, b, ptcloud_name, ptcloud_file, reconstructed_control_points, prune_flag);
%         end
% 
%         %Create a new cell array for the plyStruct2 property arrays, which contains
%         %only the reconstructed voxel x, y, z coordinates, and not the normals or
%         %colour data
%         plyStruct2.propArrayListList = cell(1, 1);
%         plyStruct2.propArrayListList{1}{1} = reconstructed_vox_pos(:, 1);   %Reconstructed voxel X coordinates
%         plyStruct2.propArrayListList{1}{2} = reconstructed_vox_pos(:, 2);   %Reconstructed voxel Y coordinates
%         plyStruct2.propArrayListList{1}{3} = reconstructed_vox_pos(:, 3);   %Reconstructed voxel Z coordinates
%         %Create a new cell array for the plyStruct2 property types, which contains 
%         %only the data types of the reconstructed voxel x, y, z coordinates, and 
%         %not the data types of the normals or colour data 
%         plyStruct2.propTypeListList = cell(1, 1);
%         plyStruct2.propTypeListList{1}(1) = "float";
%         plyStruct2.propTypeListList{1}(2) = "float";
%         plyStruct2.propTypeListList{1}(3) = "float";
%         %Create a new cell array for the plyStruct2 property names, which contains 
%         %only the names for the reconstructed voxel x, y, z coordinates, and 
%         %not for the normals or colour data 
%         plyStruct2.propNameListList = cell(1, 1);
%         plyStruct2.propNameListList{1}(1) = "x";
%         plyStruct2.propNameListList{1}(2) = "y";
%         plyStruct2.propNameListList{1}(3) = "z";
% 
%         %Write to PLY file the reconstructed voxel positions obtained from the
%         %decoder for the current input point cloud
% %         if start_lvl_cntr <= 9
% %             plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(start_lvl_cntr) '.ply'], format);
% %             disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(start_lvl_cntr) '.ply']);
% %         else
% %             plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(start_lvl_cntr) '.ply'], format);
% %             disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(start_lvl_cntr) '.ply']);
% %         end
% %         disp(' ');
%         if max_lvl_cntr <= 9
%             plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(max_lvl_cntr) '.ply'], format);
%             disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted0' num2str(max_lvl_cntr) '.ply']);
%         else
%             plyWrite(plyStruct2, [output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(max_lvl_cntr) '.ply'], format);
%             disp(['Finished writing to PLY file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_distorted' num2str(max_lvl_cntr) '.ply']);
%         end
%         disp(' ');
% 
%         %Write the geometry bitrates to the text file (the highest bitrate
%         %should come first, corresponding to the best reconstruction in
%         %..._distorted01.ply)
%         fprintf(fid_geom_bits, '%f %f\r\n', [total_geom_bpv total_geom_bits]);
%         disp(['Finished writing to geometry bitrates text file: ' output_dir ptcloud_name '_voxelized' num2str(b) '_geom_bitrates.txt']);
%         disp('------------------------------------------------------------');
% 
%         %Close all open figures
%         close all;
%         
%         %Update the counter for the next max level
%         %max_lvl_cntr = max_lvl_cntr + 1;
%     end
%    
%     %Update the counter for the next start level
%     %start_lvl_cntr = start_lvl_cntr + 1;
% end

%Close all the bitrates text files
fclose(fid_geom_bits);
if colour_compression == 1
    fclose(fid_col_bits);
    fclose(fid_all_bits);
end

%Close any other files that might still be open (for good measure)
fclose all;
    
%Finish writing to log
diary off;

disp(' ');
disp('======================= END ================================');
disp(' ');




    










