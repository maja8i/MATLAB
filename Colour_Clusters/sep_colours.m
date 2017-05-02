
%Read in input point cloud
ptCloud = pcread('swimtwoBCF.ply');
%Display input point cloud
figure;
pcshow(ptCloud, 'MarkerSize', 40);

%Choose the number of different colour clusters that we wish to have
k = 2;
%Perform k-means clustering on the colour matrix of the input point cloud
%[cluster_indices, cluster_centroids, d] = kmeans2(double(ptCloud.Color), k);
[cluster_indices, cluster_centroids, d] = kmeans(double(ptCloud.Color), k, 'Distance', 'cosine');
%Initialize a cell array to store the different sub-point-clouds resulting
%form the different colour clusters found above
sub_clouds = cell(k, 1);  
%Collect all the points belonging to the same clusters
for i = 1:k
    %Find all point row indices that belong to cluster i
    pt_rows = find(cluster_indices == i);
    %Extract the corresponding points in these rows and store them and
    %their colour values (original, for now) in the corresponding sub-cloud
    sub_clouds{i} = pointCloud(ptCloud.Location(pt_rows, :));
    sub_clouds{i}.Color = ptCloud.Color(pt_rows, :);
end

%Check total count of the points in all of the sub-clouds, to make sure 
%that they don't exceed the total number of points in the input cloud
ptcount = 0;
for j = 1:k
    ptcount = ptcount + sub_clouds{j}.Count;
end

%Plot each of the k sub-clouds separately
for pc = 1:k
    figure;
    pcshow(sub_clouds{pc}, 'MarkerSize', 40);
    %set(gca,'Color',[0 0.2 0.8]);
end

%pcwrite(sub_clouds{1}, 'BuzzPants.ply', 'PLYFormat', 'ascii');

% %Extract all unique rows from the colours matrix of ptCloud: this will tell
% %us how many unique colours are in the input point cloud
% unique_colours = unique(ptCloud.Color, 'rows'); %Colour values will be sorted in ascending order
% 
% %For each unique colour found above, find the point locations of the points
% %that have that colour and store it in the sub-cloud matrix corresponding
% %to that colour
% sub_clouds = cell(size(unique_colours, 1), 1);  %Cell array of sub-clouds divided by colour
% for i = 1:size(unique_colours, 1)
%     current_col = unique_colours(i, :);
%     %Find the corresponding points with that colour
%     inds = ismember(ptCloud.Color, current_col, 'rows');
%     pt_rows = find(inds == 1);
%     %Store the point locations and their corresponding colour in the 
%     %corresponding sub-cloud matrix, making it a pointCloud object
%     sub_clouds{i} = pointCloud(ptCloud.Location(pt_rows, :));
%     sub_clouds{i}.Color = repmat(current_col, length(pt_rows), 1);
% end