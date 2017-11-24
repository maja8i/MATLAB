%Script to compare the rate-distortion results of different codecs on one 
%or more input point clouds.

%------------------------------ User Inputs ------------------------------%

%Set the rel_deb flag to indicate whether you wish to run the 
%quality_measure C++ executable in Relase mode (rel_deb = 1) or Debug mode 
%(rel_deb = 0). Note: Release mode runs faster.  
rel_deb = 1;

%Type in the path to the quality_measure executable that you wish to use
%(use backslashes, but don't include a backslash at the end of the path 
%and don't include the executable name itself in this path)
exe_path = 'C:\Quality\quality_measure\build\Release';

%Enter the path to where the input point cloud files are contained (use
%backslashes, including a backslash at the end of the path)
data_path = '\\Pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\';

%Enter the path to the directory that contains the codec outputs (use
%backslashes, including a backslash at the end of the path)
codec_results_path = '\\Pandora\builds\test\Data\Compression\PLY\Codec_Results\';

%Enter the path to the directory where you wish to save the rate-distortion
%results (use backslashes, including a backslash at the end of the path)
rd_path = '\\Pandora\builds\test\Data\Compression\PLY\R-D_Comparisons\';

%Populate the cell array below with the names of codecs you wish to compare
codec_names = {'BezierVolume_thresh0', 'tri7', 'tri8', 'trivar'};

%Populate the cell array below with the names of input point clouds that 
%you wish to test. Don't include the _voxelizedN or .ply file extension in
%the name.
ptcloud_names = {'boxer'}; %Must be in PLY format for now

%Enter the voxelization level that you wish to use for the input (and thus
%output) point clouds. If no voxelization, put 0 here. Write the voxelizedN
%number as a string, e.g., '10'.
voxelizedN = '10';

%Populate the cell array below with numbers indicating the number of 
%reconstructed point clouds that you wish to compute distortions for, for
%each codec. This should match the number of ptcloud_name_voxelizedN_distorted#.ply 
%files in the corresponding codec_results_path\pt_cloud_name\voxelizedN\codec_name 
%directory. Note that each number must be written as a string (i.e., inside 
%quotation marks '').
nbr_reconstructions = {'5', '1', '1', '17'};    %Each column represents a different codec, each row a different input point cloud
nbr_reconstructions = repmat(nbr_reconstructions, numel(ptcloud_names), 1);    %Since each column must represent a different codec, and each row a different input point cloud (assume that a given codec produces the same no. of reconstructions for any given input point cloud)

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
all_bitrates = cell(numel(ptcloud_names), numel(codec_names)); 

%If we are using voxelized input data, compute the maximum possible
%geometric PSNR for this data (which would occur in the lossless geometry 
%reconstruction case)
if str2double(voxelizedN) > 0
    peak_val = 2^(str2double(voxelizedN)) - 1;
    max_geom_PSNR = 10*log10((peak_val^2)/(1/12));
end

%---------------------------- Run executable -----------------------------%

%For each input point cloud ...
for p = 1:numel(ptcloud_names)
    %For each codec ...
    for c = 1:numel(codec_names)
        %Obtain error measurements for this codec's reconstructions of the
        %current point cloud
        if rel_deb == 1 
            system(['cd ' exe_path ' &quality_measure.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c} ' ' codec_results_path ' ' voxelizedN]);   %Release mode
        elseif rel_deb == 0 
            system(['cd ' exe_path ' &quality_measure_d.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c} ' ' codec_results_path ' ' voxelizedN]);   %Debug mode
        end
    end %End codec loop
end %End point cloud loop

%--------------------------- Plot R-D Results ----------------------------%

%For each input point cloud ...
for p = 1:numel(ptcloud_names)
    %For each codec ...
    for c = 1:numel(codec_names)
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
        all_bitrates{p, c} = [];
        
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
        %Total (geometry + colour) bitrates
        fid = fopen([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_all_bitrates.txt']);
        if fid ~= -1
            all_bitrates{p, c} = importdata([codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_all_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' codec_results_path ptcloud_names{p} '\' vox_novox_dir '\' codec_names{c} '\' reconstruction_prefix '_all_bitrates.txt']);
        end        

        %If none of the bitrates files were able to be opened (i.e., if all of the
        %bitrates matrices above are empty), then print an error message and exit
        %the program
        if (isempty(Y_bitrates{p, c}) && isempty(U_bitrates{p, c}) && isempty(V_bitrates{p, c}) && isempty(col_bitrates{p, c}) && isempty(geom_bitrates{p, c}) && isempty(all_bitrates{p, c}))
            error('ERROR: All of the bitrates matrices for the current codec are empty. Exiting program ...');
        end
        
    end %End codec loop
    
    %Create a directory where the rate-distortion plots will be saved
    rd_dir_name = [];
    for cod = 1:numel(codec_names)
        new_addition = codec_names{cod};
        if cod == 1
            rd_dir_name = new_addition;
        else
            rd_dir_name = [rd_dir_name '_vs_' new_addition];
        end
    end
    [mkdir_success, ~, ~] = mkdir([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name]);
    if mkdir_success == 1
        disp('-------------------------------------------------------------------');
        disp(['Successfully created directory ' rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name]);
        disp('-------------------------------------------------------------------');
    else
        disp('-------------------------------------------------------------------');
        error(['ERROR: Did not create directory ' rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '. Exiting program ...']);
    end
    
    %Plot rate-distortion curves comparing the different codecs for the
    %current input point cloud
    
    %R-D plots for symmetric geometric Hausdorff distance
    figure;
    subplot(1, 2, 1);  
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 1), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);   
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 1), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_dH_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_dH_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_dH_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_dH_plots_comparison'], '-dpdf');
    
    %R-D plots for symmetric geometric RMSE
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 2), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 2), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_RMSE_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_RMSE_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_RMSE_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_RMSE_plots_comparison'], '-dpdf');

    %R-D plots for geometric PSNR
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 3), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on; 
    if (~isempty(geom_bitrates{p, i}))&&(str2double(voxelizedN) > 0)
        hold on;
        %Plot a dashed red line at the maximum possible geometric PSNR for
        %this input dataset
        %plot(bitrate_data(:, 1), repmat(max_geom_PSNR, 1, length(bitrate_data(:, 1))), 'r--');
        plot(xlim(gca), [max_geom_PSNR max_geom_PSNR], 'r--');
        legend([codec_legend_list, 'Max. possible geometric PSNR'], 'Interpreter', 'none', 'Location', 'best');
    else
        legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    end
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)'); 
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(geom_bitrates{p, i})
            bitrate_data = geom_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 3), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    if (~isempty(geom_bitrates{p, i}))&&(str2double(voxelizedN) > 0)
        hold on;
        %Plot a dashed red line at the maximum possible geometric PSNR for
        %this input dataset
        %plot(bitrate_data(:, 2), repmat(max_geom_PSNR, 1, length(bitrate_data(:, 2))), 'r--');
        plot(xlim(gca), [max_geom_PSNR max_geom_PSNR], 'r--');
        legend([codec_legend_list, 'Max. possible geometric PSNR'], 'Interpreter', 'none', 'Location', 'best');
    else
        legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    end
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_PSNR_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_PSNR_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_PSNR_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_PSNR_plots_comparison'], '-dpdf');
    
    %R-D plots for PSNR of colour Y component
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(Y_bitrates{p, i})
            bitrate_data = Y_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only Y bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(Y_bitrates{p, i})
            bitrate_data = Y_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (only Y bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_Y_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_Y_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_plots_comparison'], '-dpdf');  
    
    %R-D plots for PSNR of colour U component
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(U_bitrates{p, i})
            bitrate_data = U_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only U bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(U_bitrates{p, i})
            bitrate_data = U_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (only U bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_U_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_U_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_plots_comparison'], '-dpdf');     
    
    %R-D plots for PSNR of colour V component
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(V_bitrates{p, i})
            bitrate_data = V_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(V_bitrates{p, i})
            bitrate_data = V_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (only V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_V_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_V_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_plots_comparison'], '-dpdf');      
    
    %R-D plots for PSNR of colour Y component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_Y_vs_YUVbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_vs_YUVbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');      
    
    %R-D plots for PSNR of colour U component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_U_vs_YUVbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_vs_YUVbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');          
    
    %R-D plots for PSNR of colour V component vs (Y + U + V) bits
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(col_bitrates{p, i})
            bitrate_data = col_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_V_vs_YUVbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_vs_YUVbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');          
    
    %R-D plots for symmetric geometric Hausdorff distance vs total bitrates
    figure;
    subplot(1, 2, 1);  
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 1), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);   
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 1), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_dH_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_dH_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_dH_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_dH_vs_totalbits_plots_comparison'], '-dpdf');
    
    %R-D plots for symmetric geometric RMSE vs total bitrates
    figure;
    subplot(1, 2, 1);  
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 2), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);   
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 2), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_RMSE_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_RMSE_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_RMSE_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_RMSE_vs_totalbits_plots_comparison'], '-dpdf');
    
    %R-D plots for geometric PSNR vs total bitrates
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 3), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    if (~isempty(geom_bitrates{p, i}))&&(str2double(voxelizedN) > 0)
        hold on;
        %Plot a dashed red line at the maximum possible geometric PSNR for
        %this input dataset
        plot(bitrate_data(:, 1), repmat(max_geom_PSNR, 1, length(bitrate_data(:, 1))), 'r--');
        legend([codec_legend_list, 'Max. possible geometric PSNR'], 'Interpreter', 'none', 'Location', 'best');
    else
        legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    end
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)'); 
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 3), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    if (~isempty(geom_bitrates{p, i}))&&(str2double(voxelizedN) > 0)
        hold on;
        %Plot a dashed red line at the maximum possible geometric PSNR for
        %this input dataset
        plot(bitrate_data(:, 2), repmat(max_geom_PSNR, 1, length(bitrate_data(:, 2))), 'r--');
        legend([codec_legend_list, 'Max. possible geometric PSNR'], 'Interpreter', 'none', 'Location', 'best');
    else
        legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    end
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_geom_PSNR_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_geom_PSNR_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_PSNR_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_geom_PSNR_vs_totalbits_plots_comparison'], '-dpdf');
    
    %R-D plots for PSNR of colour Y component vs total bitrates
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 4), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_Y_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_Y_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_Y_vs_totalbits_plots_comparison'], '-dpdf');  
    
    %R-D plots for PSNR of colour U component vs total bitrates
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 5), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_U_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_U_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_U_vs_totalbits_plots_comparison'], '-dpdf');     
    
    %R-D plots for PSNR of colour V component vs total bitrates
    figure;
    subplot(1, 2, 1);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 1), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    codec_legend_list = {}; %List that will be used for the plot legend
    cll_cntr = 1;   %Counter for legend list of codec names
    for i = 1:numel(codec_names)
        if ~isempty(all_bitrates{p, i})
            bitrate_data = all_bitrates{p, i};
            errors_data = errors{p, i};
            plot(bitrate_data(:, 2), errors_data(:, 6), '-o');
            hold on;
            %Add the current codec name to the list that will be used for
            %the plot legend
            codec_legend_list{cll_cntr} = codec_names{i};
            cll_cntr = cll_cntr + 1;
        end
    end
    grid on;
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title({'Rate-Distortion Plots for', reconstruction_prefix}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_PSNR_V_vs_totalbits_plots_comparison']);
    print('-bestfit', [reconstruction_prefix '_PSNR_V_vs_totalbits_plots_comparison'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_vs_totalbits_plots_comparison']);
    print('-bestfit', [rd_path ptcloud_names{p} '\' vox_novox_dir '\' rd_dir_name '\' reconstruction_prefix '_PSNR_V_vs_totalbits_plots_comparison'], '-dpdf');    
    
    %Close all open figures before starting R-D computation for the next
    %point cloud
    close all;
    
end %End point cloud loop



