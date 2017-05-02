%Read in input point cloud
ptCloud = pcread('babyBCF.ply');

%Apply PCA on colour matrix of input point cloud (input to pca() must be a
%floating-point array, hence the conversion to single of ptCloud.Color)
[COEFF, SCORE, LATENT, ~, EXPLAINED] = pca(single(ptCloud.Color));