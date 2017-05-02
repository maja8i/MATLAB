%-------------------------------------------------------------------------%
%This function invokes a C++ program that measures geometric errors       % 
%between point clouds using several different quality metrics. Then it    %
%uses these error measurements to plot rate-distortion curves using       %
%bitrates read in from another file.                                      %
%                                                                         %          
%Inputs: rel_deb - flag to indicate whether we want to run the            %
%                  quality_measure C++ executable in Relase mode          %
%                  (rel_deb = 1) or Debug mode (rel_deb = 0)              %
%                                                                         %
%        exe_path - path to the quality_measure(_d).exe (don't include    % 
%                   the executable name in this path). Use backslashes,   %
%                   but don't include a backslash at the end of the path. %
%                                                                         %
%        data_path - path to where the point cloud files are contained.   %
%                    Use backslashes, including a backslash at the end of % 
%                    the path.                                            %
%                                                                         %
%        codec_name - a name to use on the R-D plots, so we know which    %
%                     codec the plots are for                             %
%                                                                         %
%        ptcloud_name - name of input point cloud (without the file       %
%                       extension, but assume PLY format for now) - this  %
%                       point cloud file must exist in the directory      %
%                       defined by exe_path                               %
%                                                                         %
%        nbr_reconstructions - number of reconstructed point clouds to    %
%                              compute distortions for (should match the  %
%                              number of ptcloud_name_distorted#.ply files%
%                              in the data_path\codec_name directory).    %
%                              nbr_reconstructions must be written as a   %
%                              string (i.e., in quotation marks ''), not  %
%                              as an integer (see example below).         %
%                                                                         %
%Example usage:                                                           %
%  plot_rd(1, 'C:\my_executable_releasemode_path', 'my_datapath\',        % 
%                                   'my_codec', 'babyBCF', '7')           % 
%                                                                         %
%NOTE: Currently this program expects the bitrates text files to be in the%
%      same directory as the corresponding errors text file.              %
%-------------------------------------------------------------------------%

function plot_rd(rel_deb, exe_path, data_path, codec_name, ptcloud_name, nbr_reconstructions)

%Run the quality_measure executable
if rel_deb == 1 
    system(['cd ' exe_path ' &quality_measure.exe ' data_path ' ' codec_name ' ' ptcloud_name ' ' nbr_reconstructions]);   %Release mode
elseif rel_deb == 0 
    system(['cd ' exe_path ' &quality_measure_d.exe ' data_path ' ' codec_name ' ' ptcloud_name ' ' nbr_reconstructions]); %Debug mode
end

%Open the file containing error measurements for the current codec 
errors = importdata([data_path codec_name '\' ptcloud_name '_errors.txt']);
%If "errors" is empty, then print an error message and exit the program
if isempty(errors)
    error('ERROR: "errors" is empty. Exiting program ...');
end

%Open the file(s) containing bitrates 
Y_bitrates = [];
U_bitrates = [];
V_bitrates = [];
col_bitrates = [];
geom_bitrates = [];

%Y bitrates
fid = fopen([data_path codec_name '\' ptcloud_name '_Y_bitrates.txt']);
if fid ~= -1
    Y_bitrates = importdata([data_path codec_name '\' ptcloud_name '_Y_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' data_path codec_name '\' ptcloud_name '_Y_bitrates.txt']);
end
%U bitrates
fid = fopen([data_path codec_name '\' ptcloud_name '_U_bitrates.txt']);
if fid ~= -1
    U_bitrates = importdata([data_path codec_name '\' ptcloud_name '_U_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' data_path codec_name '\' ptcloud_name '_U_bitrates.txt']);
end
%V bitrates
fid = fopen([data_path codec_name '\' ptcloud_name '_V_bitrates.txt']);
if fid ~= -1
    V_bitrates = importdata([data_path codec_name '\' ptcloud_name '_V_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' data_path codec_name '\' ptcloud_name '_V_bitrates.txt']);
end
%Total colour bitrates
fid = fopen([data_path codec_name '\' ptcloud_name '_col_bitrates.txt']);
if fid ~= -1
    col_bitrates = importdata([data_path codec_name '\' ptcloud_name '_col_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' data_path codec_name '\' ptcloud_name '_col_bitrates.txt']);
end
%Total geometry bitrates
fid = fopen([data_path codec_name '\' ptcloud_name '_geom_bitrates.txt']);
if fid ~= -1
    geom_bitrates = importdata([data_path codec_name '\' ptcloud_name '_geom_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' data_path codec_name '\' ptcloud_name '_geom_bitrates.txt']);
end

%If none of the bitrates files were able to be opened (i.e., if all of the
%bitrates matrices above are empty), then print an error message and exit
%the program
if (isempty(Y_bitrates) && isempty(U_bitrates) && isempty(V_bitrates) && isempty(col_bitrates) && isempty(geom_bitrates))
    error('ERROR: All of the bitrates matrices are empty. Exiting program ...');
end

%Plot rate-distortion curves. Here it is assumed that the errors text file
%(and therefore the "errors" matrix, above) contains, in each row, values 
%for the following error metrics, in this order: 
%   symm_hausdorff symm_rms psnr_db psnr_y psnr_u psnr_v
%Each column of "errors" is assumed to contain the error value for the
%corresponding metric for that column, for a different distorted point
%cloud (reconstruction). It is also assumed that the bitrates text files 
%(and thus each of the bitrates matrices) contain two columns of bitrates, 
%the first column containing bits per point and the second column total 
%bits in that frame, where each row of bitrates corresponds to a different 
%reconstructed version of the reference point cloud that was used to 
%compute the error values in the errors text files. Plot separate R-D 
%curves for each of the above quality metrics, against their corresponding
%bitrates.

%Geometry R-D plots
if ~isempty(geom_bitrates)
    %R-D plots for symmetric geometric Hausdorff distance
    figure;
    subplot(1, 2, 1);
    plot(geom_bitrates(:, 1), errors(:, 1), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 1), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_geom_dH_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_geom_dH_plots'], '-dpdf');

    %R-D plots for symmetric geometric RMSE
    figure;
    subplot(1, 2, 1);
    plot(geom_bitrates(:, 1), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_geom_RMSE_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_geom_RMSE_plots'], '-dpdf');

    %R-D plots for geometric PSNR
    figure;
    subplot(1, 2, 1);
    plot(geom_bitrates(:, 1), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_geom_PSNR_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_geom_PSNR_plots'], '-dpdf');
end

%Colour R-D plots
if ~isempty(Y_bitrates)
    %R-D plots for PSNR of colour Y component
    figure;
    subplot(1, 2, 1);
    plot(Y_bitrates(:, 1), errors(:, 4), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point (only Y bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    plot(Y_bitrates(:, 2), errors(:, 4), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits (only Y bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_PSNR_Y_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_PSNR_Y_plots'], '-dpdf');
end
if ~isempty(U_bitrates)
    %R-D plots for PSNR of colour U component
    figure;
    subplot(1, 2, 1);
    plot(U_bitrates(:, 1), errors(:, 5), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point (only U bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    plot(U_bitrates(:, 2), errors(:, 5), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits (only U bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_PSNR_U_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_PSNR_U_plots'], '-dpdf');
end
if ~isempty(V_bitrates)
    %R-D plots for PSNR of colour V component
    figure;
    subplot(1, 2, 1);
    plot(V_bitrates(:, 1), errors(:, 6), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point (only V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    plot(V_bitrates(:, 2), errors(:, 6), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits (only V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_PSNR_V_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_PSNR_V_plots'], '-dpdf');
end
if ~isempty(col_bitrates)
    figure;
    subplot(1, 2, 1);
    plot(col_bitrates(:, 1), errors(:, 4), '-o');
    hold on;
    plot(col_bitrates(:, 1), errors(:, 5), '-o');
    hold on;
    plot(col_bitrates(:, 1), errors(:, 6), '-o');
    grid on;
    legend('PSNR Y', 'PSNR U', 'PSNR V', 'Location', 'best');
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Bits Per Point (Y + U + V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR (dB)');
    subplot(1, 2, 2);
    plot(col_bitrates(:, 2), errors(:, 4), '-o');
    hold on;
    plot(col_bitrates(:, 2), errors(:, 5), '-o');
    hold on;
    plot(col_bitrates(:, 2), errors(:, 6), '-o');
    grid on;
    legend('PSNR Y', 'PSNR U', 'PSNR V', 'Location', 'best');
    title(['Rate-Distortion Plot for ' ptcloud_name ' using ' codec_name], 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([ptcloud_name '_' codec_name '_col_PSNR_plots']);
    print('-bestfit', [ptcloud_name '_' codec_name '_col_PSNR_plots'], '-dpdf');
end















