%-------------------------------------------------------------------------%
%This function invokes a C++ program that measures geometric errors       % 
%between point clouds using several different quality metrics. Then it    %
%uses these error measurements to plot rate-distortion curves using       %
%bitrates read in from another file.                                      %
%                                                                         %          
%----INPUTS----                                                           %
%                                                                         %
%rel_deb - flag to indicate whether we want to run the quality_measure    %
%          C++ executable in Relase mode (rel_deb = 1) or Debug mode      %
%          (rel_deb = 0)                                                  %
%                                                                         %
%exe_path - path to the quality_measure(_d).exe. Don't include the        %
%           executable name in this path. Use backslashes, but don't      %
%           include a backslash at the end of the path.                   %
%                                                                         %
%data_path - path to where the input point cloud files are. Use           %
%            backslashes, including a backslash at the end of the path.   %                                            
%                                                                         %
%codec_name - a name to use on the R-D plots, so we know which codec the  %
%             plots are for                                               %
%                                                                         %
%ptcloud_name - name of input point cloud (assume PLY format for now).    %
%               Don't include the _voxelizedN or .ply file extension.     %
%                                                                         %
%nbr_reconstructions - number of reconstructed point clouds to compute    %
%                      distortions for (should match the number of        %
%                      ptcloud_name_voxelizedN_distorted#.ply files in    % 
%                      the corresponding                                  %
%                  codec_results_path\pt_cloud_name\voxelizedN\codec_name % 
%                      directory). Note: nbr_reconstructions must be      %
%                      written as a string (i.e., in quotation marks ''), %
%                      not as an integer (see example below).             %
%                                                                         %
%codec_results_path - path to the directory that contains the codec       %
%                     outputs. Use backslashes, including a backslash at  %
%                     the end of the path.                                %
%                                                                         %
%voxelizedN - the voxelization level that you wish to use for the input   %
%             (and thus output) point clouds. If no voxelization, put 0   %
%             here. Write the voxelizedN number as a string, e.g., '10'.  %
%                                                                         %
%rd_path - path to the directory where you wish to save the R-D results   %
%          results. Use backslashes, including a backslash at the end of  %
%          the path.                                                      %
%                                                                         %
%---- EXAMPLE USAGE ----                                                  %
%                                                                         %
%  plot_rd(1, 'C:\my_executable_releasemode_path', 'my_datapath\',        % 
%                                   'my_codec', 'luna', '17',             %
%                                   'my_codec_results_path\', '10',       %
%                                   'my_rd_path\')                        % 
%                                                                         %
%-------------------------------------------------------------------------%

function plot_rd(rel_deb, exe_path, data_path, codec_name, ptcloud_name, nbr_reconstructions, codec_results_path, voxelizedN, rd_path)

%Run the quality_measure executable
if rel_deb == 1 
    system(['cd ' exe_path ' &quality_measure.exe ' data_path ' ' codec_name ' ' ptcloud_name ' ' nbr_reconstructions ' ' codec_results_path ' ' voxelizedN]);   %Release mode
elseif rel_deb == 0 
    system(['cd ' exe_path ' &quality_measure_d.exe ' data_path ' ' codec_name ' ' ptcloud_name ' ' nbr_reconstructions ' ' codec_results_path ' ' voxelizedN]);   %Release mode
end

%Open the file containing error measurements for the current codec 
if str2double(voxelizedN) == 0
    vox_novox_dir = 'OriginalPLY';
    reconstruction_prefix = ptcloud_name;
else
    vox_novox_dir = ['voxelized' voxelizedN];
    reconstruction_prefix = [ptcloud_name '_voxelized' voxelizedN];
end
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_errors.txt']);
if fid ~= -1
    errors = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_errors.txt']);
    fclose(fid);
else
    %Exit the program
    error(['ERROR: Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_errors.txt']);
end

%Open the file(s) containing bitrates 
Y_bitrates = [];
U_bitrates = [];
V_bitrates = [];
col_bitrates = [];
geom_bitrates = [];
all_bitrates = [];

%Y bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_Y_bitrates.txt']);
if fid ~= -1
    Y_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_Y_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_Y_bitrates.txt']);
end
%U bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_U_bitrates.txt']);
if fid ~= -1
    U_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_U_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_U_bitrates.txt']);
end
%V bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_V_bitrates.txt']);
if fid ~= -1
    V_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_V_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_V_bitrates.txt']);
end
%Total colour bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_col_bitrates.txt']);
if fid ~= -1
    col_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_col_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_col_bitrates.txt']);
end
%Total geometry bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_geom_bitrates.txt']);
if fid ~= -1
    geom_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_geom_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_geom_bitrates.txt']);
end
%Total (geometry + colour) bitrates
fid = fopen([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_all_bitrates.txt']);
if fid ~= -1
    all_bitrates = importdata([codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_all_bitrates.txt']);
    fclose(fid);
elseif fid == -1
    disp(['Cannot open ' codec_results_path ptcloud_name '\' vox_novox_dir '\' codec_name '\' reconstruction_prefix '_all_bitrates.txt']);
end 

%If none of the bitrates files were able to be opened (i.e., if all of the
%bitrates matrices above are empty), then print an error message and exit
%the program
if (isempty(Y_bitrates) && isempty(U_bitrates) && isempty(V_bitrates) && isempty(col_bitrates) && isempty(geom_bitrates) && isempty(all_bitrates))
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
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry only)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 1), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry only)');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_dH_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_dH_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_dH_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_dH_plots'], '-dpdf');   

    %R-D plots for symmetric geometric RMSE
    figure;
    subplot(1, 2, 1);
    plot(geom_bitrates(:, 1), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry only)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry only)');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_RMSE_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_RMSE_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_RMSE_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_RMSE_plots'], '-dpdf');   
    
    %R-D plots for geometric PSNR
    figure;
    subplot(1, 2, 1);
    plot(geom_bitrates(:, 1), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry only)');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)');
    subplot(1, 2, 2);
    plot(geom_bitrates(:, 2), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry only)');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_PSNR_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_PSNR_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_PSNR_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_PSNR_plots'], '-dpdf');      
end

%Colour R-D plots
if ~isempty(Y_bitrates)
    %R-D plots for PSNR of colour Y component
    figure;
    subplot(1, 2, 1);
    plot(Y_bitrates(:, 1), errors(:, 4), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only Y bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR Y (dB)');
    subplot(1, 2, 2);
    plot(Y_bitrates(:, 2), errors(:, 4), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (only Y bits)');   %Total bits in frame
    ylabel('Colour PSNR Y (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_PSNR_Y_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_PSNR_Y_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_Y_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_Y_plots'], '-dpdf');       
end
if ~isempty(U_bitrates)
    %R-D plots for PSNR of colour U component
    figure;
    subplot(1, 2, 1);
    plot(U_bitrates(:, 1), errors(:, 5), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only U bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR U (dB)');
    subplot(1, 2, 2);
    plot(U_bitrates(:, 2), errors(:, 5), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (only U bits)');   %Total bits in frame
    ylabel('Colour PSNR U (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_PSNR_U_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_PSNR_U_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_U_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_U_plots'], '-dpdf');           
end
if ~isempty(V_bitrates)
    %R-D plots for PSNR of colour V component
    figure;
    subplot(1, 2, 1);
    plot(V_bitrates(:, 1), errors(:, 6), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (only V bits)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR V (dB)');
    subplot(1, 2, 2);
    plot(V_bitrates(:, 2), errors(:, 6), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (only V bits)');   %Total bits in frame
    ylabel('Colour PSNR V (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_PSNR_V_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_PSNR_V_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_V_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_PSNR_V_plots'], '-dpdf');           
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
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
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
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Y + U + V bits)');   %Total bits in frame
    ylabel('Colour PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_col_PSNR_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_col_PSNR_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_col_PSNR_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_col_PSNR_plots'], '-dpdf');           
end

if ~isempty(all_bitrates)
    %R-D plots for symmetric geometric Hausdorff distance
    figure;
    subplot(1, 2, 1);
    plot(all_bitrates(:, 1), errors(:, 1), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric Hausdorff Distance');
    subplot(1, 2, 2);
    plot(all_bitrates(:, 2), errors(:, 1), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Symmetric Geometric Hausdorff Distance');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_dH_vs_totalbits_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_dH_vs_totalbits_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_dH_vs_totalbits_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_dH_vs_totalbits_plots'], '-dpdf');   

    %R-D plots for symmetric geometric RMSE
    figure;
    subplot(1, 2, 1);
    plot(all_bitrates(:, 1), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Symmetric Geometric RMSE');
    subplot(1, 2, 2);
    plot(all_bitrates(:, 2), errors(:, 2), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Symmetric Geometric RMSE');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_RMSE_vs_totalbits_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_RMSE_vs_totalbits_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_RMSE_vs_totalbits_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_RMSE_vs_totalbits_plots'], '-dpdf');   
    
    %R-D plots for geometric PSNR
    figure;
    subplot(1, 2, 1);
    plot(all_bitrates(:, 1), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Geometric PSNR (dB)');
    subplot(1, 2, 2);
    plot(all_bitrates(:, 2), errors(:, 3), '-o', 'MarkerFaceColor', 'b');
    grid on;
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Geometric PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_geom_PSNR_vs_totalbits_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_geom_PSNR_vs_totalbits_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_PSNR_vs_totalbits_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_geom_PSNR_vs_totalbits_plots'], '-dpdf'); 

    %R-D plots for colour (Y, U, V)
    figure;
    subplot(1, 2, 1);
    plot(all_bitrates(:, 1), errors(:, 4), '-o');
    hold on;
    plot(all_bitrates(:, 1), errors(:, 5), '-o');
    hold on;
    plot(all_bitrates(:, 1), errors(:, 6), '-o');
    grid on;
    legend('PSNR Y', 'PSNR U', 'PSNR V', 'Location', 'best');
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Bits Per Point (Geometry + Colour)');   %Bits per total no. of points in reconstruction
    ylabel('Colour PSNR (dB)');
    subplot(1, 2, 2);
    plot(all_bitrates(:, 2), errors(:, 4), '-o');
    hold on;
    plot(all_bitrates(:, 2), errors(:, 5), '-o');
    hold on;
    plot(all_bitrates(:, 2), errors(:, 6), '-o');
    grid on;
    legend('PSNR Y', 'PSNR U', 'PSNR V', 'Location', 'best');
    title({['Rate-Distortion Plot for ' reconstruction_prefix], ['using ' codec_name]}, 'Interpreter', 'none');
    xlabel('Total Bits (Geometry + Colour)');   %Total bits in frame
    ylabel('Colour PSNR (dB)');
    %Save the plots as a MATLAB figure and as a PDF image in the current 
    %MATLAB directory (NB: The '-bestfit' option maximizes the size of the 
    %figure to fill the page, but preserves the aspect ratio of the figure. 
    %The figure might not fill the entire page. This option leaves a 
    %minimum page margin of .25 inches).
    savefig([reconstruction_prefix '_' codec_name '_col_PSNR_vs_totalbits_plots']);
    print('-bestfit', [reconstruction_prefix '_' codec_name '_col_PSNR_vs_totalbits_plots'], '-dpdf');
    %Save in our network directory as well
    savefig([rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_col_PSNR_vs_totalbits_plots']);
    print('-bestfit', [rd_path ptcloud_name '\' vox_novox_dir '\' reconstruction_prefix '_' codec_name '_col_PSNR_vs_totalbits_plots'], '-dpdf');           
end


