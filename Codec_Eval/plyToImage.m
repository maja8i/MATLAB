function cdata = plyToImage(infile, degrees, N)
%PLYTOIMAGE creates a NxNx3 image from ply file
%
% infile = ply file to read
% degrees (optional) = degrees to rotate ply file around
%   clockwise looking down the Y axis [default 0]
% N (optional) = image width [default width of PLY file, else 1024]
%
% cdata = NxNx3 image

% Process optional args.
if nargin < 3
    N = 1024;
    if nargin < 2
        degrees = 0;
    end
end

% Construct point cloud from file.
vpc = voxelizedPointCloud(infile);

% Convert from world to frame coordinates if necessary.
if vpc.FrameToWorldScale == 1
    % In world coordinates
    cubeWidth = N-1;
    vpc = vpc.setTransform(cubeWidth);
    vpc = vpc.worldToFrame();
end

% Create Morton Codes and sort.
vpc = vpc.mortonizeAndSort();

% Voxelize.
vpc = vpc.voxelize();

% Snapshot.
cdata = vpc.orthogonalProjection(degrees);

end

