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
data_path = 'C:\Quality\quality_measure\data\';

%Populate the cell array below with the names of codecs you wish to compare
codec_names = {'RAHT_USQ_RLGR', 'HVR53', 'HEVC'};

%Populate the cell array below with the names of input point clouds
%(without the file type extension) that you wish to test
ptcloud_names = {'luna_voxelized10'}; %Must be in PLY format for now

%Populate the cell array below with numbers indicating the number of 
%reconstructed point clouds that you wish to compute distortions for, for
%each codec. This should match the number of ptcloud_name_distorted#.ply 
%files in the corresponding data_path\codec_name directory. Note that each
%number must be written as a string (i.e., inside quotation marks '').
nbr_reconstructions = {'17', '11', '11'};    %Each column represents a different codec, each row a different input point cloud

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

%For each input point cloud ...
for p = 1:numel(ptcloud_names)
    %For each codec ...
    for c = 1:numel(codec_names)
        %Obtain error measurements for this codec's reconstructions of the
        %current point cloud
        if rel_deb == 1 
            system(['cd ' exe_path ' &quality_measure.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c}]);   %Release mode
        elseif rel_deb == 0 
            system(['cd ' exe_path ' &quality_measure_d.exe ' data_path ' ' codec_names{c} ' ' ptcloud_names{p} ' ' nbr_reconstructions{p, c}]); %Debug mode
        end
        
        %Open the file containing error measurements for the current codec 
        errors{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_errors.txt']);
        %If "errors" is empty, then print an error message and exit the program
        if isempty(errors{p, c})
            error('ERROR: "errors" is empty. Exiting program ...');
        end

        %Open the file(s) containing bitrates 
        Y_bitrates{p, c} = [];
        U_bitrates{p, c} = [];
        V_bitrates{p, c} = [];
        col_bitrates{p, c} = [];
        geom_bitrates{p, c} = [];
        
        %Y bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_Y_bitrates.txt']);
        if fid ~= -1
            Y_bitrates{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_Y_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_Y_bitrates.txt']);
        end
        %U bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_U_bitrates.txt']);
        if fid ~= -1
            U_bitrates{p, c} = importdata([data_path '\' codec_names{c} '\' ptcloud_names{p} '_U_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_U_bitrates.txt']);
        end
        %V bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_V_bitrates.txt']);
        if fid ~= -1
            V_bitrates{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_V_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_V_bitrates.txt']);
        end
        %Total colour bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
        if fid ~= -1
            col_bitrates{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_col_bitrates.txt']);
        end
        %Total geometry bitrates
        fid = fopen([data_path codec_names{c} '\' ptcloud_names{p} '_geom_bitrates.txt']);
        if fid ~= -1
            geom_bitrates{p, c} = importdata([data_path codec_names{c} '\' ptcloud_names{p} '_geom_bitrates.txt']);
            fclose(fid);
        elseif fid == -1
            disp(['Cannot open ' data_path codec_names{c} '\' ptcloud_names{p} '_geom_bitrates.txt']);
        end

        %If none of the bitrates files were able to be opened (i.e., if all of the
        %bitrates matrices above are empty), then print an error message and exit
        %the program
        if (isempty(Y_bitrates{p, c}) && isempty(U_bitrates{p, c}) && isempty(V_bitrates{p, c}) && isempty(col_bitrates{p, c}) && isempty(geom_bitrates{p, c}))
            error('ERROR: All of the bitrates matrices for the current codec are empty. Exiting program ...');
        end
        
    end %End codec loop
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_geom_dH_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_geom_dH_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_geom_RMSE_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_geom_RMSE_plots_comparison'], '-dpdf');
    
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
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    legend(codec_legend_list, 'Interpreter', 'none', 'Location', 'best');
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_geom_PSNR_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_geom_PSNR_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (only Y bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_Y_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_Y_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (only U bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_U_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_U_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (only V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_V_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_V_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_Y_vs_YUVbits_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_Y_vs_YUVbits_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_U_vs_YUVbits_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_U_vs_YUVbits_plots_comparison'], '-dpdf');
    
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
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
    title(['Rate-Distortion Plots for ' ptcloud_names{p}], 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_names{p} '_PSNR_V_vs_YUVbits_plots_comparison']);
    print('-bestfit', [ptcloud_names{p} '_PSNR_V_vs_YUVbits_plots_comparison'], '-dpdf');
    
    %Close all open figures before starting R-D computation for the next
    %point cloud
    close all;
    
end %End point cloud loop



