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

%-------------------------------------------------------------------------%

%Initialize the errors and bitrates cell arrays, to hold the error values
%and bitrates for the different point clouds and codecs. For each cell
%array below, each point cloud is a new cell row and each codec is a new 
%cell column. 
errors = cell(numel(ptcloud_names), numel(codec_names));   
col_bitrates = cell(numel(ptcloud_names), numel(codec_names));

%For each input point cloud
for p = 1:numel(ptcloud_names)
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
    
    %For each codec
    for c = 1:numel(codec_names)
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
        lbri_cntr = lbri_cntr + 1;
        
        %In the col_bitrates text file for the current point cloud and
        %codec, find the closest bpp value to "high_bitrate" (i.e., find
        %its row number in the col_bitrates text file)
        [diff2, highb_index] = min(abs(col_bitrates{p, c}(:, 1) - high_bitrate));
        high_bitrate_reconstruction_inds(hbri_cntr) = highb_index;
        hbri_cntr = hbri_cntr + 1;
        
     end     %End codec loop
     
     
    %Initialize vectors to store axis handles for the different subplots
    %that will be generated below
    h_axes_low = zeros(1, numel(codec_names));   %For low bitrate subplots
    h_axes_high = zeros(1, numel(codec_names));   %For high bitrate subplots
        
    %Each element (column) of low_bitrate_reconstruction_inds contains
    %the index of the closest bitrate to low_bitrate for a different 
    %codec. These indices correspond to reconstruction numbers of the 
    %reconstructed PLY files. Pull out these corresponding PLY 
    %reconstructions, compute orthographic projections of them from the 
    %same viewpoint, and plot these projections on the same figure. 
    figure;
    %Plot original PLY
    h_axes_low(1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), 1);
    [plyStruct, A] = plyRead([data_path ptcloud_names{p} '.ply']); %Read in PLY file
    col_data = [A(:, 4), A(:, 5), A(:, 6)]; %Extract PLY colour data
    scatter3(A(:, 1), A(:, 2), A(:, 3), 10, col_data./255);    %Plot point cloud (divide col_data elements by 255 because RGB values in MATLAB range from 0 to 1)
    axis equal tight;   %Maintain aspect ratio of point cloud
    grid off;
    axis off;
    ax = gca;
    ax.Clipping = 'off';   %Turn off axis clipping (do don't get holes when zoom in)
    title(['Original ' ptcloud_names{p} '.ply'], 'Interpreter', 'none');
    %Plot each of the different codec reconstructions
    for i = 1:length(low_bitrate_reconstruction_inds)
        h_axes_low(i+1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), i+1);
        if (low_bitrate_reconstruction_inds(i) <= 9)
            [plyStruct, A] = plyRead([data_path codec_names{i} '\' ptcloud_names{p} '_distorted0' num2str(i) '.ply']); %Read in PLY file
        else
            [plyStruct, A] = plyRead([data_path codec_names{i} '\' ptcloud_names{p} '_distorted' num2str(i) '.ply']); %Read in PLY file
        end
        col_data = [A(:, 4), A(:, 5), A(:, 6)]; %Extract PLY colour data
        scatter3(A(:, 1), A(:, 2), A(:, 3), 10, col_data./255);    %Plot point cloud (divide col_data elements by 255 because RGB values in MATLAB range from 0 to 1)
        axis equal tight;   %Maintain aspect ratio of point cloud
        grid off;
        axis off;
        ax = gca;
        ax.Clipping = 'off';   %Turn off axis clipping (do don't get holes when zoom in)
        title([codec_names{i} ' Reconstruction at ' num2str(col_bitrates{p, i}(low_bitrate_reconstruction_inds(i), 1)) ' bpp (Y+U+V)'], 'Interpreter', 'none');
    end
    %Synchronize axes properties of the different subplots
    linkprop(h_axes_low, {'CameraPosition', 'CameraUpVector', 'View'});
    
    %Each element (column) of high_bitrate_reconstruction_inds contains
    %the index of the closest bitrate to high_bitrate for a different 
    %codec. These indices correspond to reconstruction numbers of the 
    %reconstructed PLY files. Pull out these corresponding PLY 
    %reconstructions, compute orthographic projections of them from the 
    %same viewpoint, and plot these projections on the same figure.  
    figure;
    %Plot original PLY
    h_axes_high(1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), 1);
    [plyStruct, A] = plyRead([data_path ptcloud_names{p} '.ply']); %Read in PLY file
    col_data = [A(:, 4), A(:, 5), A(:, 6)]; %Extract PLY colour data
    scatter3(A(:, 1), A(:, 2), A(:, 3), 10, col_data./255);    %Plot point cloud (divide col_data elements by 255 because RGB values in MATLAB range from 0 to 1)
    axis equal tight;   %Maintain aspect ratio of point cloud
    grid off;
    axis off;
    ax = gca;
    ax.Clipping = 'off';   %Turn off axis clipping (do don't get holes when zoom in)
    title(['Original ' ptcloud_names{p} '.ply'], 'Interpreter', 'none');
    %Plot each of the different codec reconstructions
    for i = 1:length(high_bitrate_reconstruction_inds)
        h_axes_high(i+1) = subplot(floor((numel(codec_names) + 1)/2), ceil((numel(codec_names) + 1)/2), i+1);
        if (high_bitrate_reconstruction_inds(i) <= 9)
            [plyStruct, A] = plyRead([data_path codec_names{i} '\' ptcloud_names{p} '_distorted0' num2str(i) '.ply']); %Read in PLY file
        else
            [plyStruct, A] = plyRead([data_path codec_names{i} '\' ptcloud_names{p} '_distorted' num2str(i) '.ply']); %Read in PLY file
        end
        col_data = [A(:, 4), A(:, 5), A(:, 6)]; %Extract PLY colour data
        scatter3(A(:, 1), A(:, 2), A(:, 3), 10, col_data./255);    %Plot point cloud (divide col_data elements by 255 because RGB values in MATLAB range from 0 to 1)
        axis equal tight;   %Maintain aspect ratio of point cloud
        grid off;
        axis off;
        ax = gca;
        ax.Clipping = 'off';   %Turn off axis clipping (do don't get holes when zoom in)
        title([codec_names{i} ' Reconstruction at ' num2str(col_bitrates{p, i}(high_bitrate_reconstruction_inds(i), 1)) ' bpp (Y+U+V)'], 'Interpreter', 'none');
    end
    %Synchronize axes properties of the different subplots
    linkprop(h_axes_low, {'CameraPosition', 'CameraUpVector', 'View'});
       
end     %End point cloud loop