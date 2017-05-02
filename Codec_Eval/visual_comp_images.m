%Run this code from inside the "Codec_Eval" directory.

%Requires directory "Phil": add this directory plus its sub-directories to
%the current MATLAB path.

%Also requires directory "subplot_tight": add this directory plus its 
%sub-directories to the current MATLAB path.

%------------------------------ User Inputs ------------------------------%

%Enter the path to where the input point cloud files are contained (use
%backslashes, including a backslash at the end of the path)
data_path = '\\Pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\';

%Enter the path to the directory that contains the codec outputs (use
%backslashes, including a backslash at the end of the path)
codec_results_path = '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\';

%Populate the cell array below with the names of codecs whose results you
%wish to compare
codec_names = {'HEVC_INTRA', 'CubeProj_HEVC_INTRA_RedundancyRemoval', 'CubeProj_HEVC_INTRA_SCD', 'CubeProj_HEVC_INTRA_SCD_FillingInterp', 'CubeProj_HEVC_INTRA_SCD_RedundancyRemoval'};

%Populate the cell array below with the names of input point clouds that 
%you wish to test. Don't include the _voxelizedN or .ply file extension in
%the name.
ptcloud_names = {'redandblack_dress'}; %Must be in PLY format for now

%Enter the voxelization level that you wish to use for the input (and thus
%output) point clouds. If no voxelization, put 0 here. Write the voxelizedN
%number as a string, e.g., '10';
voxelizedN = '10';

%Enter the path to the directory where you wish to save the rate-distortion
%results (use backslashes, including a backslash at the end of the path)
rd_path = '\\Pandora\builds\test\Data\Compression\PLY\R-D_Comparisons\';

%Enter the name of the directory inside rd_path, where you wish to store
%the images of the reconstructions and the SSIM images that will be
%obtained from this M-file. NOTE: IT IS ASSUMED THAT THIS DIRECTORY ALREADY
%EXISTS.
%rd_dir_name = 'RAHT_USQ_RLGR_vs_HEVC_INTRA_vs_Mathematica_v02';

%Enter a low bitrate (bits per point) that you wish to compare 
%reconstructions for
low_bitrate = 1.0;
%Enter a high bitrate (bits per point) that you wish to compare 
%reconstructions for
high_bitrate = 3.5;

%Choose degrees that you wish to rotate PLY files around to get orthogonal
%image projections
degrees = [0, 90, 180, 270];

%-------------------------------------------------------------------------%

%Get the directory name inside rd_path, where the images of the 
%reconstructions and the SSIM images will be stored. **NOTE: IT IS ASSUMED 
%THAT THIS DIRECTORY ALREADY EXISTS, SO IT IS NOT CREATED HERE.** 
rd_dir_name = [];
for cod = 1:numel(codec_names)
    new_addition = codec_names{cod};
    if cod == 1
        rd_dir_name = new_addition;
    else
        rd_dir_name = [rd_dir_name '_vs_' new_addition];
    end
end
%rd_dir_name = 'CubeProj_HEVC_INTRA';

%Initialize the errors and bitrates cell arrays, to hold the error values
%and bitrates for the different point clouds and codecs. For each cell
%array below, each point cloud is a new cell row and each codec is a new 
%cell column. 
errors = cell(numel(ptcloud_names), numel(codec_names));   
col_bitrates = cell(numel(ptcloud_names), numel(codec_names));

%For each input point cloud
for p = 1:numel(ptcloud_names)
    disp('==============================================================');
    disp(['Processing ' ptcloud_names{p} '.ply ...']);
    %Initialize an array to store the indices of the reconstructed
    %point clouds that we wish to compare for all the codecs,
    %corresponding to the current input point cloud and "low_bitrate"
    %reconstructions
    low_bitrate_reconstruction_inds = [];
    lbri_cntr = 1;
    %Initialize an array to store the indices of the reconstructed
    %point clouds that we wish to compare for all the codecs,
    %corresponding to the current input point cloud and "high_bitrate"
    %reconstructions
    high_bitrate_reconstruction_inds = [];
    hbri_cntr = 1;
    %Initialize cell arrays to store images of the current input point
    %cloud
    input_imgs = cell(length(degrees), 1);  %RGB space
    input_imgs_yuv = cell(length(degrees), 1);  %YUV space
    input_imgs_y2rgb_disp = cell(length(degrees), 1);   %Y image converted to RGB with U and V set to 0, for display purposes only
    input_imgs_u2rgb_disp = cell(length(degrees), 1);   %U image converted to RGB with Y and V set to 0, for display purposes only
    input_imgs_v2rgb_disp = cell(length(degrees), 1);   %V image converted to RGB with U and Y set to 0, for display purposes only
    %Initialize cell arrays to store images of reconstructions at a low
    %bitrate, for each codec
    recon_imgs_low = cell(length(degrees), numel(codec_names)); %RGB space
    recon_imgs_low_yuv = cell(length(degrees), numel(codec_names)); %YUV space
    recon_imgs_low_y2rgb_disp = cell(length(degrees), numel(codec_names));   %Y image converted to RGB with U and V set to 0, for display purposes only
    recon_imgs_low_u2rgb_disp = cell(length(degrees), numel(codec_names));   %U image converted to RGB with Y and V set to 0, for display purposes only
    recon_imgs_low_v2rgb_disp = cell(length(degrees), numel(codec_names));   %V image converted to RGB with U and Y set to 0, for display purposes only    
    %Initialize cell arrays to store images of reconstructions at a high
    %bitrate, for each codec
    recon_imgs_high = cell(length(degrees), numel(codec_names));    %RGB space
    recon_imgs_high_yuv = cell(length(degrees), numel(codec_names));    %YUV space
    recon_imgs_high_y2rgb_disp = cell(length(degrees), numel(codec_names));   %Y image converted to RGB with U and V set to 0, for display purposes only
    recon_imgs_high_u2rgb_disp = cell(length(degrees), numel(codec_names));   %U image converted to RGB with Y and V set to 0, for display purposes only
    recon_imgs_high_v2rgb_disp = cell(length(degrees), numel(codec_names));   %V image converted to RGB with U and Y set to 0, for display purposes only        
    %Initialize arrays to store SSIM values and difference images (local
    %SSIM maps) for images of reconstructions at a low bitrate
    ssim_low_vals = [];
    ssim_low_maps = cell(length(degrees), numel(codec_names));
    %Initialize arrays to store SSIM values and difference images (local
    %SSIM maps) for images of reconstructions at a high bitrate  
    ssim_high_vals = [];
    ssim_high_maps = cell(length(degrees), numel(codec_names));
    %Initialize arrays to store PSNR values for reconstructed images at a
    %low bitrate
    psnr_low_imgs_y = [];
    psnr_low_imgs_u = [];
    psnr_low_imgs_v = [];
    %Initialize arrays to store PSNR values for reconstructed images at a
    %high bitrate
    psnr_high_imgs_y = [];
    psnr_high_imgs_u = [];
    psnr_high_imgs_v = [];
    %Initialize arrays to store average PSNR values for low and high
    %bitrates for each codec for the current point cloud 
    avg_psnr_low_y = [];
    avg_psnr_low_u = [];
    avg_psnr_low_v = [];
    avg_psnr_high_y = [];
    avg_psnr_high_u = [];
    avg_psnr_high_v = [];
    %Initialize arrays to store average SSIM values for low and high
    %bitrates for each codec for the current point cloud 
    avg_ssim_low = [];
    avg_ssim_high = [];
        
    %Create directory names that will be used later
    if str2double(voxelizedN) == 0
        vox_novox_dir = 'OriginalPLY';
        reconstruction_prefix = ptcloud_names{p};
    else
        vox_novox_dir = ['voxelized' voxelizedN];
        reconstruction_prefix = [ptcloud_names{p} '_voxelized' voxelizedN];
    end
    
    %Check if the directory rd_dir_name exists; if not, print an error 
    %message and exit the program
   
    exist_res = exist([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\'], 'dir');
    if exist_res == 7
        disp('-------------------------------------------------------------');
        disp(['Directory ' rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name ' exists']);
    elseif exist_res == 0
        disp('-------------------------------------------------------------');
        error(['ERROR: Directory ' rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name ' does not exist. Exiting program ...']);
    end
    
    %Get images of input point cloud at different viewpoints
    for i = 1:length(degrees)
        disp('-------------------------------------------------------------');
        disp(['Rendering input image at degree ' num2str(degrees(i)) ' ...']);
        input_imgs{i, 1} = plyToImage([data_path vox_novox_dir '\' reconstruction_prefix '.ply'], degrees(i));
    end
    
    %For each codec
    for c = 1:numel(codec_names)
        disp('-------------------------------------------------------------');
        disp(['Processing ' codec_names{c} ' ...']);
        %Open the file containing error measurements for the current codec 
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
        if fid ~= -1
            errors{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
            fclose(fid);
        else
            %Exit the program
            error(['ERROR: Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
        end

        %Open the file containing bitrates 
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
        if fid ~= -1
            col_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
        end
        
        %In the col_bitrates text file for the current point cloud and
        %codec, find the closest bpp value to "low_bitrate" (i.e., find
        %its row number in the col_bitrates text file)
        [diff1, lowb_index] = min(abs(col_bitrates{p, c}(:, 1) - low_bitrate));
        low_bitrate_reconstruction_inds(lbri_cntr) = lowb_index;
        %Get images of the reconstruction corresponding to the above low
        %bitrate, for the current codec, at different viewpoints
        for i = 1:length(degrees)
            disp('-------------------------------------------------------------');
            disp(['Rendering image for low bitrate at degree ' num2str(degrees(i)) ' ...']);
            if (low_bitrate_reconstruction_inds(lbri_cntr) <= 9)
                recon_imgs_low{i, c} = plyToImage([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_distorted0' num2str(low_bitrate_reconstruction_inds(lbri_cntr)) '.ply'], degrees(i));
            else
                recon_imgs_low{i, c} = plyToImage([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_distorted' num2str(low_bitrate_reconstruction_inds(lbri_cntr)) '.ply'], degrees(i));
            end
        end
        %Update lbri_cntr
        lbri_cntr = lbri_cntr + 1;
        
        %In the col_bitrates text file for the current point cloud and
        %codec, find the closest bpp value to "high_bitrate" (i.e., find
        %its row number in the col_bitrates text file)
        [diff2, highb_index] = min(abs(col_bitrates{p, c}(:, 1) - high_bitrate));
        high_bitrate_reconstruction_inds(hbri_cntr) = highb_index;
        %Get images of the reconstruction corresponding to the above high
        %bitrate, for the current codec, at different viewpoints
        for i = 1:length(degrees)
            disp('-------------------------------------------------------------');
            disp(['Rendering image for high bitrate at degree ' num2str(degrees(i)) ' ...']);
            if (high_bitrate_reconstruction_inds(hbri_cntr) <= 9)
                recon_imgs_high{i, c} = plyToImage([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_distorted0' num2str(high_bitrate_reconstruction_inds(hbri_cntr)) '.ply'], degrees(i));
            else
                recon_imgs_high{i, c} = plyToImage([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_distorted' num2str(high_bitrate_reconstruction_inds(hbri_cntr)) '.ply'], degrees(i));
            end
        end
        %Update hbri_cntr
        hbri_cntr = hbri_cntr + 1;
   
     end     %End codec loop
     
     %----------------------- RGB to YUV Conversion ----------------------%
     
     for i = 1:length(degrees)
         %Convert the input image at the current viewpoint to YUV space
         input_imgs_yuv{i} = rgbToYuvImage(input_imgs{i}); 
         %Convert each codec reconstruction image at the corresponding
         %viewpoint to YUV space
         for j = 1:numel(codec_names)
             %For low bitrate
             recon_imgs_low_yuv{i, j} = rgbToYuvImage(recon_imgs_low{i, j}); 
             %For high bitrate
             recon_imgs_high_yuv{i, j} = rgbToYuvImage(recon_imgs_high{i, j});
         end
     end
     
%      %----------------------- YUV to RGB Conversion ----------------------%
%      
%      %This is for display purposes only: so that MATLAB can display the Y,
%      %U, and V images in a conventional way
%      
%      for i = 1:length(degrees)
%          %Input Y image conversion
%          temp_input = input_imgs_yuv{i};
%          temp_input(:, :, 2:3) = zeros; %Set U and V components to zeros
%          input_imgs_y2rgb_disp{i} = yuvToRgbImage(temp_input);
%          %Convert to uint8 so MATLAB can display the image
%          input_imgs_y2rgb_disp{i} = uint8(input_imgs_y2rgb_disp{i});
%          %Input U image conversion
%          temp_input = input_imgs_yuv{i};
%          temp_input(:, :, 1) = zeros; %Set Y component to zeros
%          temp_input(:, :, 3) = zeros; %Set V component to zeros
%          input_imgs_u2rgb_disp{i} = yuvToRgbImage(temp_input);   
%          %Convert to uint8 so MATLAB can display the image
%          input_imgs_u2rgb_disp{i} = uint8(input_imgs_u2rgb_disp{i});         
%          %Input V image conversion
%          temp_input = input_imgs_yuv{i};
%          temp_input(:, :, 1:2) = zeros; %Set Y and U components to zeros
%          input_imgs_v2rgb_disp{i} = yuvToRgbImage(temp_input);  
%          %Convert to uint8 so MATLAB can display the image
%          input_imgs_v2rgb_disp{i} = uint8(input_imgs_v2rgb_disp{i});         
%          
%          %Convert each Y/U/V codec reconstruction image at the current
%          %viewpoint to RGB space
%          for j = 1:numel(codec_names)
%              %For low bitrate ...
%              
%              %Y image conversion
%              temp_input = recon_imgs_low_yuv{i, j};
%              temp_input(:, :, 2:3) = zeros; %Set U and V components to zeros
%              recon_imgs_low_y2rgb_disp{i, j} = yuvToRgbImage(temp_input);
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_low_y2rgb_disp{i, j} = uint8(recon_imgs_low_y2rgb_disp{i, j});
%              %U image conversion
%              temp_input = recon_imgs_low_yuv{i, j};
%              temp_input(:, :, 1) = zeros; %Set Y component to zeros
%              temp_input(:, :, 3) = zeros; %Set V component to zeros
%              recon_imgs_low_u2rgb_disp{i, j} = yuvToRgbImage(temp_input);
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_low_u2rgb_disp{i, j} = uint8(recon_imgs_low_u2rgb_disp{i, j});
%              %V image conversion
%              temp_input = recon_imgs_low_yuv{i, j};
%              temp_input(:, :, 1:2) = zeros; %Set Y and U components to zeros
%              recon_imgs_low_v2rgb_disp{i, j} = yuvToRgbImage(temp_input);
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_low_v2rgb_disp{i, j} = uint8(recon_imgs_low_v2rgb_disp{i, j});
%           
%              %For high bitrate ...
% 
%              %Y image conversion
%              temp_input = recon_imgs_high_yuv{i, j};
%              temp_input(:, :, 2:3) = zeros; %Set U and V components to zeros
%              recon_imgs_high_y2rgb_disp{i, j} = yuvToRgbImage(temp_input);
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_high_y2rgb_disp{i, j} = uint8(recon_imgs_high_y2rgb_disp{i, j});         
%              %U image conversion
%              temp_input = recon_imgs_high_yuv{i, j};
%              temp_input(:, :, 1) = zeros; %Set Y component to zeros
%              temp_input(:, :, 3) = zeros; %Set V component to zeros
%              recon_imgs_high_u2rgb_disp{i, j} = yuvToRgbImage(temp_input);
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_high_u2rgb_disp{i, j} = uint8(recon_imgs_high_u2rgb_disp{i, j});
%              %V image conversion
%              temp_input = recon_imgs_high_yuv{i, j};
%              temp_input(:, :, 1:2) = zeros; %Set Y and U components to zeros
%              recon_imgs_high_v2rgb_disp{i, j} = yuvToRgbImage(temp_input);  
%              %Convert to uint8 so MATLAB can display the image
%              recon_imgs_high_v2rgb_disp{i, j} = uint8(recon_imgs_high_v2rgb_disp{i, j});
%          end
%      end
     
     %---------------------------- Low bitrate ---------------------------%
     
     %Display all the image reconstructions obtained from the different 
     %codecs at a low bitrate, side by side, in RGB space
     for i = 1:length(degrees)
         figure;
         %Display input image at the current degree value
         ax_low(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs{i});
         axis equal;
         axis off;
         title(['Original ' reconstruction_prefix '.ply'], 'Interpreter', 'none');
         %Display images of each of the different codec reconstructions at 
         %the current degree value 
         for j = 1:numel(codec_names)
            ax_low(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_low{i, j});
            axis equal;
            axis off;  
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)']}, 'Interpreter', 'none');
         end
         linkaxes(ax_low);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [reconstruction_prefix '_images_low_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' reconstruction_prefix '_images_low_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');      
     end
     
     %Display all the Y components of the image reconstructions obtained 
     %from the different codecs at a low bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input Y image at the current degree value
         ax_low_y(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 1), [0 255]); colormap('gray');
         %imagesc(input_imgs_y2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'Y Component'}, 'Interpreter', 'none');
         %Display Y images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_low_y(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_low_yuv{i, j}(:, :, 1), [0 255]); colormap('gray');
            %imagesc(recon_imgs_low_y2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute Y PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_low_imgs_y(i, j) = psnr(uint8(recon_imgs_low_yuv{i, j}(:, :, 1)), uint8(input_imgs_yuv{i}(:, :, 1)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR Y = ' num2str(psnr_low_imgs_y(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_low_y);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_Y_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_Y_low_bitrate_' num2str(degrees(i)) 'degrees']);
     end
   
     %Display all the U components of the image reconstructions obtained 
     %from the different codecs at a low bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input U image at the current degree value
         ax_low_u(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 2), [0 255]); colormap('winter');
         %imagesc(input_imgs_u2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'U Component'}, 'Interpreter', 'none');
         %Display U images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_low_u(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_low_yuv{i, j}(:, :, 2), [0 255]); colormap('winter');
            %imagesc(recon_imgs_low_u2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute U PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_low_imgs_u(i, j) = psnr(uint8(recon_imgs_low_yuv{i, j}(:, :, 2)), uint8(input_imgs_yuv{i}(:, :, 2)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR U = ' num2str(psnr_low_imgs_u(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_low_u);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_U_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_U_low_bitrate_' num2str(degrees(i)) 'degrees']);
     end
     
     %Display all the V components of the image reconstructions obtained 
     %from the different codecs at a low bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input V image at the current degree value
         ax_low_v(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 3), [0 255]); colormap('autumn');
         %imagesc(input_imgs_v2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'V Component'}, 'Interpreter', 'none');
         %Display V images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_low_v(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_low_yuv{i, j}(:, :, 3), [0 255]); colormap('autumn');
            %imagesc(recon_imgs_low_v2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute V PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_low_imgs_v(i, j) = psnr(uint8(recon_imgs_low_yuv{i, j}(:, :, 3)), uint8(input_imgs_yuv{i}(:, :, 3)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR V = ' num2str(psnr_low_imgs_v(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_low_v);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_V_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_V_low_bitrate_' num2str(degrees(i)) 'degrees']);
     end
     
     %Compute SSIM between the input point cloud images and each of the
     %corresponding reconstructed images from the different codecs, in RGB
     %space
     disp('-------------------------------------------------------------');
     disp('Computing SSIM for low-bitrate images ...');
     for i = 1:length(degrees)
         figure;
         for j = 1:numel(codec_names)
            %Compute SSIM between the current reconstruction image and the
            %corresponding original image
            [ssim_low_vals(i, j), ssim_low_maps{i, j}] = ssim(recon_imgs_low{i, j}, input_imgs{i});
            %Display SSIM difference images for all the different codec
            %reconstructions side by side
            ax_ssim_low(j) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j);
            imshow(ssim_low_maps{i, j});
            title({[codec_names{j} ' at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)'], ['Global (mean) SSIM = ' num2str(ssim_low_vals(i, j))]}, 'Interpreter', 'none');
         end
         linkaxes(ax_ssim_low);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_SSIM_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [reconstruction_prefix '_SSIM_low_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_SSIM_low_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' reconstruction_prefix '_SSIM_low_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');      
     end
     
     %--------------------------- High bitrate ---------------------------%
     
     %Display all the image reconstructions obtained from the different 
     %codecs at a high bitrate, side by side, in RGB space
     for i = 1:length(degrees)
         figure;
         %Display input image at the current degree value
         ax_high(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs{i});
         axis equal;
         axis off;
         title(['Original ' reconstruction_prefix '.ply'], 'Interpreter', 'none');
         %Display images of each of the different codec reconstructions at 
         %the current degree value 
         for j = 1:numel(codec_names)
            ax_high(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_high{i, j});
            axis equal;
            axis off;
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)']}, 'Interpreter', 'none');
         end
         linkaxes(ax_high);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_high_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [reconstruction_prefix '_images_high_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_high_bitrate_' num2str(degrees(i)) 'degrees']);
         %print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' reconstruction_prefix '_images_high_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');      
     end
     
     %Display all the Y components of the image reconstructions obtained 
     %from the different codecs at a high bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input Y image at the current degree value
         ax_high_y(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 1), [0 255]); colormap('gray');
         %imagesc(input_imgs_y2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'Y Component'}, 'Interpreter', 'none');
         %Display Y images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_high_y(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_high_yuv{i, j}(:, :, 1), [0 255]); colormap('gray');
            %imagesc(recon_imgs_high_y2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute Y PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_high_imgs_y(i, j) = psnr(uint8(recon_imgs_high_yuv{i, j}(:, :, 1)), uint8(input_imgs_yuv{i}(:, :, 1)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR Y = ' num2str(psnr_high_imgs_y(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_high_y);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_Y_high_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_Y_high_bitrate_' num2str(degrees(i)) 'degrees']);
     end
   
     %Display all the U components of the image reconstructions obtained 
     %from the different codecs at a high bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input U image at the current degree value
         ax_high_u(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 2), [0 255]); colormap('winter');
         %imagesc(input_imgs_u2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'U Component'}, 'Interpreter', 'none');
         %Display U images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_high_u(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_high_yuv{i, j}(:, :, 2), [0 255]); colormap('winter');
            %imagesc(recon_imgs_high_u2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute U PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_high_imgs_u(i, j) = psnr(uint8(recon_imgs_high_yuv{i, j}(:, :, 2)), uint8(input_imgs_yuv{i}(:, :, 2)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR U = ' num2str(psnr_high_imgs_u(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_high_u);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_U_high_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_U_high_bitrate_' num2str(degrees(i)) 'degrees']);
     end
     
     %Display all the V components of the image reconstructions obtained 
     %from the different codecs at a high bitrate, side by side
     for i = 1:length(degrees)         
         figure;
         %Display input V image at the current degree value
         ax_high_v(1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), 1);
         imagesc(input_imgs_yuv{i}(:, :, 3), [0 255]); colormap('autumn');
         %imagesc(input_imgs_v2rgb_disp{i});
         axis equal;
         axis off;
         title({['Original ' reconstruction_prefix '.ply'], 'V Component'}, 'Interpreter', 'none');
         %Display V images of each of the different codec reconstructions
         %at the current degree value 
         for j = 1:numel(codec_names)
            ax_high_v(j+1) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j+1);
            imagesc(recon_imgs_high_yuv{i, j}(:, :, 3), [0 255]); colormap('autumn');
            %imagesc(recon_imgs_high_v2rgb_disp{i, j});
            axis equal;
            axis off;
            %Compute V PSNR between the current reconstruction image and
            %the corresponding original image
            psnr_high_imgs_v(i, j) = psnr(uint8(recon_imgs_high_yuv{i, j}(:, :, 3)), uint8(input_imgs_yuv{i}(:, :, 3)));
            title({[codec_names{j} ' Reconstruction'], ['at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V); PSNR V = ' num2str(psnr_high_imgs_v(i, j)) ' dB']}, 'Interpreter', 'none');
         end
         linkaxes(ax_high_v);
         %Save the images in the current MATLAB directory
         savefig([reconstruction_prefix '_images_V_high_bitrate_' num2str(degrees(i)) 'degrees']);
         %Save the images in our network directory as well
         savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_images_V_high_bitrate_' num2str(degrees(i)) 'degrees']);
     end     
     
     %Compute SSIM between the input point cloud images and each of the
     %corresponding reconstructed images from the different codecs, in RGB
     %space
     disp('-------------------------------------------------------------');
     disp('Computing SSIM for high-bitrate images ...');
     for i = 1:length(degrees)
        figure;
        for j = 1:numel(codec_names)
            %Compute SSIM between the current reconstruction image and the
            %corresponding original image
            [ssim_high_vals(i, j), ssim_high_maps{i, j}] = ssim(recon_imgs_high{i, j}, input_imgs{i});                     
            %Display SSIM difference images for all the different codec
            %reconstructions side by side
            ax_ssim_high(j) = subplot_tight(ceil((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 2)/2), j);
            imshow(ssim_high_maps{i, j});
            title({[codec_names{j} ' at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)'], ['Global (mean) SSIM = ' num2str(ssim_high_vals(i, j))]}, 'Interpreter', 'none');            
        end
        linkaxes(ax_ssim_high);
        %Save the images in the current MATLAB directory
        savefig([reconstruction_prefix '_SSIM_high_bitrate_' num2str(degrees(i)) 'degrees']);
        %print('-bestfit', [reconstruction_prefix '_SSIM_high_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');
        %Save the images in our network directory as well
        savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_SSIM_high_bitrate_' num2str(degrees(i)) 'degrees']);
        %print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' reconstruction_prefix '_SSIM_high_bitrate_' num2str(degrees(i)) 'degrees'], '-dpdf');      
     end       
     
     %----------------------------- R-D Plots ----------------------------%
     
     %Compute the average PSNR Y, U, and V values from images at all 
     %viewpoints that were used, for low and high bitrates separately, for
     %each codec
     for cod = 1:numel(codec_names)
         %Low bitrate
         avg_psnr_low_y(cod) = mean(psnr_low_imgs_y(:, cod));
         avg_psnr_low_u(cod) = mean(psnr_low_imgs_u(:, cod));
         avg_psnr_low_v(cod) = mean(psnr_low_imgs_v(:, cod));
         %High bitrate
         avg_psnr_high_y(cod) = mean(psnr_high_imgs_y(:, cod));
         avg_psnr_high_u(cod) = mean(psnr_high_imgs_u(:, cod));
         avg_psnr_high_v(cod) = mean(psnr_high_imgs_v(:, cod));
     end
     PSNR_data_y = [avg_psnr_low_y; avg_psnr_high_y];
     PSNR_data_u = [avg_psnr_low_u; avg_psnr_high_u];
     PSNR_data_v = [avg_psnr_low_v; avg_psnr_high_v];
     
     %Compute the average SSIM value from images at all viewpoints that
     %were used, for low and high bitrates separately, for each codec
     for cod = 1:numel(codec_names)
         avg_ssim_low(cod) = mean(ssim_low_vals(:, cod));
         avg_ssim_high(cod) = mean(ssim_high_vals(:, cod));
     end
     SSIM_data = [avg_ssim_low; avg_ssim_high];
     
     %Image PSNR Y plots
     figure;
     codec_legend_list = {}; %List that will be used for the plot legend
     cll_cntr = 1;   %Counter for legend list of codec names
     for i = 1:numel(codec_names)
        bitrates = [col_bitrates{p, i}(low_bitrate_reconstruction_inds(i), 1); col_bitrates{p, i}(high_bitrate_reconstruction_inds(i), 1)];
        plot(bitrates, PSNR_data_y(:, i), '-o');
        hold on;
        %Add the current codec name to the list that will be used for
        %the plot legend
        codec_legend_list{cll_cntr} = codec_names{i};
        cll_cntr = cll_cntr + 1;
     end
     grid on;
     legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
     title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
     xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
     ylabel('Average Image PSNR Y (dB)');     
     %Save the plots as a MATLAB figure and as a PDF image in the current 
     %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
     %figure to fill the page, but preserves the aspect ratio of the figure. 
     %The figure might not fill the entire page. This option leaves a 
     %minimum page margin of .25 inches).
     savefig([reconstruction_prefix '_image_PSNR_Y_plots_comparison']);
     print('-bestfit', [reconstruction_prefix '_image_PSNR_Y_plots_comparison'], '-dpdf');
     %Save in our network directory as well
     savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_Y_plots_comparison']);
     print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_Y_plots_comparison'], '-dpdf');              
     
     %Image PSNR U plots
     figure;
     codec_legend_list = {}; %List that will be used for the plot legend
     cll_cntr = 1;   %Counter for legend list of codec names
     for i = 1:numel(codec_names)
        bitrates = [col_bitrates{p, i}(low_bitrate_reconstruction_inds(i), 1); col_bitrates{p, i}(high_bitrate_reconstruction_inds(i), 1)];
        plot(bitrates, PSNR_data_u(:, i), '-o');
        hold on;
        %Add the current codec name to the list that will be used for
        %the plot legend
        codec_legend_list{cll_cntr} = codec_names{i};
        cll_cntr = cll_cntr + 1;
     end
     grid on;
     legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
     title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
     xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
     ylabel('Average Image PSNR U (dB)');     
     %Save the plots as a MATLAB figure and as a PDF image in the current 
     %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
     %figure to fill the page, but preserves the aspect ratio of the figure. 
     %The figure might not fill the entire page. This option leaves a 
     %minimum page margin of .25 inches).
     savefig([reconstruction_prefix '_image_PSNR_U_plots_comparison']);
     print('-bestfit', [reconstruction_prefix '_image_PSNR_U_plots_comparison'], '-dpdf');
     %Save in our network directory as well
     savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_U_plots_comparison']);
     print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_U_plots_comparison'], '-dpdf');              

     %Image PSNR V plots
     figure;
     codec_legend_list = {}; %List that will be used for the plot legend
     cll_cntr = 1;   %Counter for legend list of codec names
     for i = 1:numel(codec_names)
        bitrates = [col_bitrates{p, i}(low_bitrate_reconstruction_inds(i), 1); col_bitrates{p, i}(high_bitrate_reconstruction_inds(i), 1)];
        plot(bitrates, PSNR_data_v(:, i), '-o');
        hold on;
        %Add the current codec name to the list that will be used for
        %the plot legend
        codec_legend_list{cll_cntr} = codec_names{i};
        cll_cntr = cll_cntr + 1;
     end
     grid on;
     legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
     title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
     xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
     ylabel('Average Image PSNR V (dB)');     
     %Save the plots as a MATLAB figure and as a PDF image in the current 
     %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
     %figure to fill the page, but preserves the aspect ratio of the figure. 
     %The figure might not fill the entire page. This option leaves a 
     %minimum page margin of .25 inches).
     savefig([reconstruction_prefix '_image_PSNR_V_plots_comparison']);
     print('-bestfit', [reconstruction_prefix '_image_PSNR_V_plots_comparison'], '-dpdf');
     %Save in our network directory as well
     savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_V_plots_comparison']);
     print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_image_PSNR_V_plots_comparison'], '-dpdf');              

     %Image SSIM plots
     figure;
     codec_legend_list = {}; %List that will be used for the plot legend
     cll_cntr = 1;   %Counter for legend list of codec names
     for i = 1:numel(codec_names)
        bitrates = [col_bitrates{p, i}(low_bitrate_reconstruction_inds(i), 1); col_bitrates{p, i}(high_bitrate_reconstruction_inds(i), 1)];
        plot(bitrates, SSIM_data(:, i), '-o');
        hold on;
        %Add the current codec name to the list that will be used for
        %the plot legend
        codec_legend_list{cll_cntr} = codec_names{i};
        cll_cntr = cll_cntr + 1;
     end
     grid on;
     legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
     title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
     xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
     ylabel('Average Image SSIM');    
     %Save the plots as a MATLAB figure and as a PDF image in the current 
     %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
     %figure to fill the page, but preserves the aspect ratio of the figure. 
     %The figure might not fill the entire page. This option leaves a 
     %minimum page margin of .25 inches).
     savefig([reconstruction_prefix '_SSIM_plots_comparison']);
     print('-bestfit', [reconstruction_prefix '_SSIM_plots_comparison'], '-dpdf');
     %Save in our network directory as well
     savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_SSIM_plots_comparison']);
     print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_SSIM_plots_comparison'], '-dpdf');          
    
     %Close all open figures before processing the next point cloud
     close all;
     
end     %End point cloud loop