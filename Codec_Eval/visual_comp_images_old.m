%------------------------------ User Inputs ------------------------------%

%Enter the path to where the point cloud files are contained (use
%backslashes, including a backslash at the end of the path). It is assumed
%that this is also the root directory for the bitrates and errors text 
%files for each codec and point cloud.
data_path = 'C:\Quality\quality_measure\data\';

%Populate the cell array below with the names of codecs whose results you
%wish to compare
codec_names = {'RAHT_USQ_RLGR', 'HVR53', 'HEVC'};

%Populate the cell array below with the names of input point clouds
%(without the file type extension) whose results you wish to compare
ptcloud_names = {'luna_voxelized10'}; %Must be in PLY format for now

%Enter a low bitrate (bits per point) that you wish to compare 
%reconstructions for
low_bitrate = 1.4;
%Enter a high bitrate (bits per point) that you wish to compare 
%reconstructions for
high_bitrate = 3.4;

%Choose degrees that you wish to rotate PLY files around to get orthogonal
%image projections
degrees = [180, 90];

%-------------------------------------------------------------------------%

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
    %Initialize a cell array to store images of the current input point
    %cloud
    input_imgs = cell(length(degrees), 1);
    %Initialize a cell array to store images of reconstructions at a low
    %bitrate, for each codec
    recon_imgs_low = cell(length(degrees), numel(codec_names));
    %Initialize a cell array to store images of reconstructions at a high
    %bitrate, for each codec
    recon_imgs_high = cell(length(degrees), numel(codec_names));
    %Initialize cell arrays to store SSIM values and difference images
    %(local SSIM maps) for images of reconstructions at a low bitrate
    ssim_low_vals = cell(length(degrees), numel(codec_names));
    ssim_low_maps = cell(length(degrees), numel(codec_names));
    %Initialize cell arrays to store SSIM values and difference images
    %(local SSIM maps) for images of reconstructions at a high bitrate  
    ssim_high_vals = cell(length(degrees), numel(codec_names));
    ssim_high_maps = cell(length(degrees), numel(codec_names));
    
    %Get images of input point cloud at different viewpoints
    for i = 1:length(degrees)
        disp('-------------------------------------------------------------');
        disp(['Rendering input image at degree ' num2str(degrees(i)) ' ...']);
        input_imgs{i, 1} = plyToImage([data_path ptcloud_names{p} '.ply'], degrees(i));
    end
    
    %For each codec
    for c = 1:numel(codec_names)
        disp('-------------------------------------------------------------');
        disp(['Processing ' codec_names{c} ' ...']);
        %Open the file containing error measurements for the current codec 
        errors{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_errors.txt']);
        %If "errors" is empty, then print an error message and exit the program
        if isempty(errors{p, c})
            error('ERROR: "errors" is empty. Exiting program ...');
        end

        %Open the file containing bitrates 
        col_bitrates{p, c} = [];    %Total colour bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
        if fid ~= -1
            col_bitrates{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            error(['ERROR: Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
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
                recon_imgs_low{i, c} = plyToImage([data_path codec_names{c} '\' ptcloud_names{p} '_distorted0' num2str(low_bitrate_reconstruction_inds(lbri_cntr)) '.ply'], degrees(i));
            else
                recon_imgs_low{i, c} = plyToImage([data_path codec_names{c} '\' ptcloud_names{p} '_distorted' num2str(low_bitrate_reconstruction_inds(lbri_cntr)) '.ply'], degrees(i));
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
                recon_imgs_high{i, c} = plyToImage([data_path codec_names{c} '\' ptcloud_names{p} '_distorted0' num2str(high_bitrate_reconstruction_inds(hbri_cntr)) '.ply'], degrees(i));
            else
                recon_imgs_high{i, c} = plyToImage([data_path codec_names{c} '\' ptcloud_names{p} '_distorted' num2str(high_bitrate_reconstruction_inds(hbri_cntr)) '.ply'], degrees(i));
            end
        end
        %Update hbri_cntr
        hbri_cntr = hbri_cntr + 1;
   
     end     %End codec loop
     
     %---------------------------- Low bitrate ---------------------------%
     
     %Display all the image reconstructions obtained from the different 
     %codecs at a low bitrate, side by side
     for i = 1:length(degrees)
         figure;
         %Display input image at the current degree value
         ax_low(1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), 1);
         image(input_imgs{i});
         axis equal;
         axis off;
         title(['Original ' ptcloud_names{p} '.ply'], 'Interpreter', 'none');
         %Display images of each of the different codec reconstructions at 
         %the current degree value 
         for j = 1:numel(codec_names)
            ax_low(j+1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), j+1);
            image(recon_imgs_low{i, j});
            axis equal;
            axis off;
            title([codec_names{j} ' Reconstruction at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)'], 'Interpreter', 'none');
         end
         linkaxes(ax_low);
     end

     %Compute SSIM between the input point cloud images and each of the
     %corresponding reconstructed images from the different codecs
     disp('-------------------------------------------------------------');
     disp('Computing SSIM for low-bitrate images ...');
     for i = 1:length(degrees)
         figure;
         for j = 1:numel(codec_names)
            [ssim_low_vals{i, j}, ssim_low_maps{i, j}] = ssim(recon_imgs_low{i, j}, input_imgs{i});
            %Display SSIM difference images for all the different codec
            %reconstructions side by side
            ax_ssim_low(j) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), j);
            imshow(ssim_low_maps{i, j});
            title({[codec_names{j} ' Reconstruction at ' num2str(col_bitrates{p, j}(low_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)']; ['Global (mean) SSIM = ' num2str(ssim_low_vals{i, j})]}, 'Interpreter', 'none');
         end
         linkaxes(ax_ssim_low);
     end

     
     %--------------------------- High bitrate ---------------------------%
     
     %Display all the image reconstructions obtained from the different 
     %codecs at a high bitrate, side by side
     for i = 1:length(degrees)
         figure;
         %Display input image at the current degree value
         ax_high(1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), 1);
         image(input_imgs{i});
         axis equal;
         axis off;
         title(['Original ' ptcloud_names{p} '.ply'], 'Interpreter', 'none');
         %Display images of each of the different codec reconstructions at 
         %the current degree value 
         for j = 1:numel(codec_names)
            ax_high(j+1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), j+1);
            image(recon_imgs_high{i, j});
            axis equal;
            axis off;
            title([codec_names{j} ' Reconstruction at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)'], 'Interpreter', 'none');
         end
         linkaxes(ax_high);
     end
     
     %Compute SSIM between the input point cloud images and each of the
     %corresponding reconstructed images from the different codecs
     disp('-------------------------------------------------------------');
     disp('Computing SSIM for high-bitrate images ...');
     for i = 1:length(degrees)
        figure;
        for j = 1:numel(codec_names)
            %Display SSIM difference images for all the different codec
            %reconstructions side by side
            [ssim_high_vals{i, j}, ssim_high_maps{i, j}] = ssim(recon_imgs_high{i, j}, input_imgs{i});
            ax_ssim_high(j) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), j);
            imshow(ssim_high_maps{i, j});
            title({[codec_names{j} ' Reconstruction at ' num2str(col_bitrates{p, j}(high_bitrate_reconstruction_inds(j), 1)) ' bpp (Y+U+V)']; ['Global (mean) SSIM = ' num2str(ssim_high_vals{i, j})]}, 'Interpreter', 'none');
        end
        linkaxes(ax_ssim_high);
     end       
end     %End point cloud loop