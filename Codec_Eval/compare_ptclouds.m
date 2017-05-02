%Script to compare the rate-distortion results for the same codec on
%different input point clouds. 

%------------------------------ User Inputs ------------------------------%

%Set the rel_deb flag to indicate whether you wish to run the 
%quality_measure C++ executable in Relase mode (rel_deb = 1) or Debug mode 
%(rel_deb = 0). Note: Release mode runs faster.  
rel_deb = 1;

%Type in the path to the quality_measure executable that you wish to use
%(use backslashes, but don't include a backslash at the end of the path 
%and don't include the executable name itself in this path)
exe_path = 'C:\Quality\quality_measure\build\Release';

%Enter the path to where the point cloud files are contained (use
%backslashes, including a backslash at the end of the path)
data_path = '\\Pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\';

%Enter the path to the directory that contains the codec outputs (use
%backslashes, including a backslash at the end of the path)
codec_results_path = '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\';

%Enter the path to the directory where you wish to save the rate-distortion
%results (use backslashes, including a backslash at the end of the path)
rd_path = '\\Pandora\builds\test\Data\Compression\PLY\R-D_Comparisons\';

%Populate the cell array below with the names of codecs you wish to get
%results for
codec_names = {'RAHT_USQ_RLGR', 'HEVC_INTRA'};

%Populate the cell array below with the names of input point clouds that 
%you wish to test. Don't include the _voxelizedN or .ply file extension in
%the name.
ptcloud_names = {'luna'}; %Must be in PLY format for now

%Enter the voxelization level that you wish to use for the input (and thus
%output) point clouds. If no voxelization, put 0 here. Write the voxelizedN
%number as a string, e.g., '10';
voxelizedN = '10';

%Populate the cell array below with numbers indicating the number of 
%reconstructed point clouds that you wish to compute distortions for, for
%each codec. This should match the number of ptcloud_name_voxelizedN_distorted#.ply 
%files in the corresponding codec_results_path\pt_cloud_name\voxelizedN\codec_name 
%directory. Note that each number must be written as a string (i.e., inside 
%quotation marks '').
%nbr_reconstructions = {'17', '11'; '17', '11'};    %Each column represents a different codec, each row a different input point cloud
nbr_reconstructions = {'17', '11'};    %Each column represents a different codec, each row a different input point cloud

%-------------------------------------------------------------------------%

%Initialize the errors and bitrates cell arrays, to hold the error values
%and bitrates for the different point clouds and codecs. For each cell
%array below, each point cloud is a new cell row and each codec is a new 
%cell column. 
errors = cell(numel(ptcloud_names), numel(codec_names));   
geom_bitrates = cell(numel(ptcloud_names), numel(codec_names)); 
Y_bitrates = cell(numel(ptcloud_names), numel(codec_names));  
U_bitrates = cell(numel(ptcloud_names), numel(codec_names)); 
V_bitrates = cell(numel(ptcloud_names), numel(codec_names)); 
col_bitrates = cell(numel(ptcloud_names), numel(codec_names)); 

%For each codec ...
for c = 1:numel(codec_names)
    %For each input point cloud ...
    for p = 1:numel(ptcloud_names)
        %Obtain error measurements for this codec's reconstructions of the
        %current point cloud
        if rel_deb == 1 
            system(['cd ' exe_path ' &quality_measure.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c} ' ' codec_results_path ' ' voxelizedN]);   %Release mode
        elseif rel_deb == 0 
            system(['cd ' exe_path ' &quality_measure_d.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c} ' ' codec_results_path ' ' voxelizedN]);   %Release mode
        end
        
        %Open the file containing error measurements for the current codec 
        if str2double(voxelizedN) == 0
            vox_novox_dir = 'OriginalPLY';
            reconstruction_prefix = ptcloud_names{p};
        else
            vox_novox_dir = ['voxelized' voxelizedN];
            reconstruction_prefix = [ptcloud_names{p} '_voxelized' voxelizedN];
        end
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
        if fid ~= -1
            errors{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
            fclose(fid);
        else
            %Exit the program
            error(['ERROR: Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_errors.txt']);
        end

        %Open the file(s) containing bitrates 
        Y_bitrates{p, c} = [];
        U_bitrates{p, c} = [];
        V_bitrates{p, c} = [];
        col_bitrates{p, c} = [];
        geom_bitrates{p, c} = [];
        
        %Y bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_Y_bitrates.txt']);
        if fid ~= -1
            Y_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_Y_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_Y_bitrates.txt']);
        end
        %U bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_U_bitrates.txt']);
        if fid ~= -1
            U_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_U_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_U_bitrates.txt']);
        end
        %V bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_V_bitrates.txt']);
        if fid ~= -1
            V_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_V_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_V_bitrates.txt']);
        end
        %Total colour bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
        if fid ~= -1
            col_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_col_bitrates.txt']);
        end
        %Total geometry bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_geom_bitrates.txt']);
        if fid ~= -1
            geom_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_geom_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_geom_bitrates.txt']);
        end

        %If none of the bitrates files were able to be opened (i.e., if all of the
        %bitrates matrices above are empty), then print an error message and exit
        %the program
        if (isempty(Y_bitrates{p, c}) && isempty(U_bitrates{p, c}) && isempty(V_bitrates{p, c}) && isempty(col_bitrates{p, c}) && isempty(geom_bitrates{p, c}))
            error('ERROR: All of the bitrates matrices for the current codec are empty. Exiting program ...');
        end
        
    end %End point cloud loop
    
    %Plot rate-distortion curves comparing the results for the same codec
    %on different input point clouds 
    
    %R-D plots for symmetric geometric Hausdorff distance
    figure;
    subplot(1, 2, 1);  
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 1), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);   
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 1), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as 
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_geom_dH_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_geom_dH_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_geom_dH_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_geom_dH_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_dH_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_dH_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_dH_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_dH_plots_comparison'], '-dpdf');    
    end
    
    %R-D plots for symmetric geometric RMSE
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 2), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if voxelizedN == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 2), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_geom_RMSE_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_geom_RMSE_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_geom_RMSE_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_geom_RMSE_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_RMSE_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_RMSE_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_RMSE_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_RMSE_plots_comparison'], '-dpdf');    
    end   
    
    %R-D plots for geometric PSNR
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 3), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)'); 
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(geom_bitrates{i, c})
            bitrate_data = geom_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 3), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_geom_PSNR_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_geom_PSNR_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_geom_PSNR_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_geom_PSNR_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_PSNR_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_PSNR_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_PSNR_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_geom_PSNR_plots_comparison'], '-dpdf');    
    end   
    
    %R-D plots for PSNR of colour Y component
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(Y_bitrates{i, c})
            bitrate_data = Y_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 4), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (only Y bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(Y_bitrates{i, c})
            bitrate_data = Y_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 4), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (only Y bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_Y_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_Y_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_Y_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_Y_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_plots_comparison'], '-dpdf');    
    end 
   
    %R-D plots for PSNR of colour U component
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(U_bitrates{i, c})
            bitrate_data = U_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 5), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (only U bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(U_bitrates{i, c})
            bitrate_data = U_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 5), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (only U bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_U_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_U_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_U_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_U_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_plots_comparison'], '-dpdf');    
    end    
    
    %R-D plots for PSNR of colour V component
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(V_bitrates{i, c})
            bitrate_data = V_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 6), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (only V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(V_bitrates{i, c})
            bitrate_data = V_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 6), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (only V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_V_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_V_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_V_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_V_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_plots_comparison'], '-dpdf');    
    end      
    
    %R-D plots for PSNR of colour Y component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 4), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 4), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_Y_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_Y_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');    
    end   
    
    %R-D plots for PSNR of colour U component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 5), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 5), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_U_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_U_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');    
    end     
    
    %R-D plots for PSNR of colour V component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 1), errors_data(:, 6), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    ptcloud_legend_list = {}; %List that will be used for the plot legend
    pll_cntr = 1;   %Counter for legend list of point cloud names
    for i = 1:numel(ptcloud_names)
        if ~isempty(col_bitrates{i, c})
            bitrate_data = col_bitrates{i, c};
            errors_data = errors{i, c};
            plot(bitrate_data(:, 2), errors_data(:, 6), '-o');
            hold on;
            %Add the current point cloud name to the list that will be used
            %for the plot legend
            ptcloud_legend_list{pll_cntr} = ptcloud_names{i};
            pll_cntr = pll_cntr + 1;
        end
    end
    grid on;
    legend(ptcloud_legend_list, 'Interpreter', 'none', 'Location', 'best');
    if str2double(voxelizedN) == 0
        title(['Rate-Distortion Plots for ' codec_names{c}], 'Interpreter', 'none');
    else
        title({['Rate-Distortion Plots for ' codec_names{c}], ['using voxelized' voxelizedN ' Versions of Input Clouds']}, 'Interpreter', 'none');
    end
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches). Save in our network directory as
    %well.
    if str2double(voxelizedN) == 0
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_PSNR_V_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_PSNR_V_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');      
    else
        %Current MATLAB directory
        savefig([codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_vs_YUVbits_plots_comparison']);
        print('-bestfit', [codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');    
        %Network directory
        savefig([rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_vs_YUVbits_plots_comparison']);
        print('-bestfit', [rd_path codec_names{c} '_diff_ptclouds_voxelized' voxelizedN '_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');    
    end     

    %Close all open figures before starting R-D computation for the next
    %codec
    close all;
    
end %End codec loop



