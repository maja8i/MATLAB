%Requires directory "MinBoundSphere&Circle": add this directory plus its 
%sub-directories to the current MATLAB path.

%Requires directory "Phil": add this directory plus its sub-directories to
%the current MATLAB path.

%Read in input PLY point cloud
[~, ptcloud, ~] = plyRead('\\pandora\builds\test\Data\Compression\PLY\Point_Clouds\8i\voxelized10\gladiator_voxelized10.ply');

%Find the minimum bounding sphere for the input point cloud
[r, C, Xb] = ExactMinBoundSphere3D(ptcloud(:, 1:3));
%Display the input point cloud in its minimum bounding sphere
[~] = VisualizeBoundSphere(ptcloud(:, 1:3), r, C, ptcloud(:, 4:6)./255);
hold on;

%Extract the x, y, and z coordinates of the centre of the bounding sphere
xc = C(1);
yc = C(2);
zc = C(3);

%Initialize a matrix to store the projected coordinates of the input point
%cloud (projected onto the surface of the bounding sphere)
proj_coords = zeros(size(ptcloud, 1), 3);

%For each input (x, y, z) point
for i = 1:size(ptcloud, 1)
    %Write the point in a coordinate system centred at the centre of the
    %bounding sphere
    p_temp(1, 1) = ptcloud(i, 1) - xc;
    p_temp(1, 2) = ptcloud(i, 2) - yc;
    p_temp(1, 3) = ptcloud(i, 3) - zc;
    %Compute the length of the p_temp vector
    norm_p_temp = norm(p_temp);
    %Scale p_temp so that its length is equal to the radius of the bounding
    %sphere
    scaled_p_temp = (r/norm_p_temp).*p_temp;
    %Compute the projection (coordinates) of the input point on the surface
    %of the bounding sphere
    proj_coords(i, 1) = scaled_p_temp(1) + xc;
    proj_coords(i, 2) = scaled_p_temp(2) + yc;
    proj_coords(i, 3) = scaled_p_temp(3) + zc;
end

%Plot the projected points on the surface of the minimum bounding box of
%the input point cloud
scatter3(proj_coords(:,1), proj_coords(:,2), proj_coords(:,3), 5, [ptcloud(:, 4)./255, ptcloud(:, 5)./255, ptcloud(:, 6)./255], 'filled');
hold off;

%Reconstruct original points from projected points



