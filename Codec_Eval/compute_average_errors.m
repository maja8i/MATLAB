%Script to compute the average error values across all frames of a given
%dynamic point cloud sequence, for a number of different distortion levels. 

%IMPORTANT: Assumes that the error values for individual frames at each
%distortion level have already been computed and written to file.

% %------------------------------ User Inputs ------------------------------%
% 
% %Enter the path to the directory that contains the codec outputs (use
% %backslashes, including a backslash at the end of the path)
% codec_results_path = '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\Dynamic\';
% 
% %Enter the name of the input point cloud. Don't include the _voxN or .ply 
% %file extension in the name
% ptcloud_name = 'redandblack'; %Must be in PLY format for now
% 
% %Enter the voxelization level for the input point cloud. If no voxelization, 
% %put 0 here.
% voxN = 10;
% 
% %Enter the codec name for which you wish to compute the average errors
% codec_name = 'RAHT_USQ_RLGR';
% 
% %Enter the number of different distortion levels whose errors you wish to
% %use in the computation of the average. This number cannot exceed the 
% %number of "distorted" folders corresponding to the selected point cloud
% %and codec. 
% dist_levels = 17;    

function [] = compute_average_errors(codec_results_path, ptcloud_name, voxN, codec_name, dist_levels)

%-------------------------------------------------------------------------%

if voxN == 0
    voxdir = '\';
    voxname = '_';
else
    voxdir = ['\vox' num2str(voxN) '\'];
    voxname = ['_vox' num2str(voxN) '_'];
end

%Initialize a matrix to hold the average error values for each distortion 
%level 
avg_errors = zeros(dist_levels, 6);    %6 average error values per distortion level
 
%For each distortion level
for d_lvl = 1:dist_levels
    if d_lvl <= 9
        distorted_foldername = [ptcloud_name voxname 'distorted0' num2str(d_lvl)];
    else
        distorted_foldername = [ptcloud_name voxname 'distorted' num2str(d_lvl)];
    end
    
    %Open the errors text file at that distortion level, which contains
    %errors for each frame
    fid = fopen([codec_results_path ptcloud_name voxdir codec_name '\' distorted_foldername '\' distorted_foldername '_errors.txt']);
    if fid ~= -1
        errors = importdata([codec_results_path ptcloud_name voxdir codec_name '\' distorted_foldername '\' distorted_foldername '_errors.txt']);
        fclose(fid);
    else
        %Exit the program
        error(['ERROR: Cannot open ' codec_results_path ptcloud_name voxdir codec_name '\' distorted_foldername '\' distorted_foldername '_errors.txt']);
    end
     
    %Compute the average of all the errors found in "errors", and store
    %them in their correct location inside avg_errors
    avg_errors(d_lvl, :) = sum(errors, 1)./size(errors, 1);
end

%Open an errors text file in the location where the avg_matrix data
%is supposed to be saved, ready for writing (this will overwrite any
%existing data in that file, if it already exists, which is what we want to
%happen)
fid2 = fopen([codec_results_path ptcloud_name voxdir codec_name '\' ptcloud_name voxname  'errors.txt'], 'w');
if fid2 ~= -1
    %Write the avg_matrix data to this file, in the same format (i.e., 6
    %error values per distortion level)
    for d_lvl = 1:size(avg_errors, 1)
        fprintf(fid2, '%f %f %f %f %f %f\r\n', avg_errors(d_lvl, 1), avg_errors(d_lvl, 2), avg_errors(d_lvl, 3), avg_errors(d_lvl, 4), avg_errors(d_lvl, 5), avg_errors(d_lvl, 6));
    end
    fclose(fid2);
else
    %Exit the program
    error(['ERROR: Cannot open ' codec_results_path ptcloud_name voxdir codec_name '\' ptcloud_name voxname  'errors.txt'])
end

