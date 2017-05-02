% Test projection

% infile = 'D:\Data\Compression\PLY\raw_inputs\longdress1180.ply';
infile = 'D:\Data\Compression\PLY\voxelized10\longdress1180_voxelized10.ply';
degrees = 90;

% Create image.
cdata = plyToImage(infile,degrees);

% Display.
image(cdata);
axis equal;
