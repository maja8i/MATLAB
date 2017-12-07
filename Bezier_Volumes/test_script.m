%Corners of leaf cell 5171 (its index before pruning; after pruning, its
%index is 3540) at level 8 (root = level 1)
cnrs = [159.5000  431.5000  159.5000
  167.5000  431.5000  159.5000
  167.5000  439.5000  159.5000
  159.5000  439.5000  159.5000
  159.5000  431.5000  167.5000
  167.5000  431.5000  167.5000
  167.5000  439.5000  167.5000
  159.5000  439.5000  167.5000];

hold on gcf;
scatter3(cnrs(:, 1), cnrs(:, 2), cnrs(:, 3), 30, 'filled', 'm');

%Reconstructed voxels associated with the above leaf cell (when don't set
%control points inside (-q/2, q/2) to 0)
vox = [164	438	160
	165	439	160
	164	439	160
	164	439	161
	160	436	160
	161	437	160
	160	437	160
	162	437	160
	162	438	160
	163	438	160
	163	439	160
	162	439	160
	162	438	161
	163	439	161
	162	439	161
	160	438	160
	161	438	160
	160	438	161
	161	438	161
	161	439	161
	160	439	161
	161	439	162
	160	439	162];

hold on gcf;
scatter3(vox(:, 1), vox(:, 2), vox(:, 3), 10, 'filled', 'k');

%Reconstructed control points
cp = [2.2266
    4.1309
    0.6787
   -1.7773
    9.0234
   10.5322
    2.0825
    4.6289];

%Write the control point value next to the corresponding corner
hold on gcf; text(cnrs(1, 1), cnrs(1, 2), cnrs(1, 3), num2str(cp(1)));
hold on gcf; text(cnrs(2, 1), cnrs(2, 2), cnrs(2, 3), num2str(cp(2)));
hold on gcf; text(cnrs(3, 1), cnrs(3, 2), cnrs(3, 3), num2str(cp(3)));
hold on gcf; text(cnrs(4, 1), cnrs(4, 2), cnrs(4, 3), num2str(cp(4)));
hold on gcf; text(cnrs(5, 1), cnrs(5, 2), cnrs(5, 3), num2str(cp(5)));
hold on gcf; text(cnrs(6, 1), cnrs(6, 2), cnrs(6, 3), num2str(cp(6)));
hold on gcf; text(cnrs(7, 1), cnrs(7, 2), cnrs(7, 3), num2str(cp(7)));
hold on gcf; text(cnrs(8, 1), cnrs(8, 2), cnrs(8, 3), num2str(cp(8)));

% %Original voxels associated with the above leaf cell
% orig = [160	435	160
% 	161	435	160
% 	162	435	160
% 	163	435	160
% 	160	436	160
% 	160	436	161
% 	160	437	161
% 	161	436	160
% 	161	436	161
% 	161	437	161
% 	160	438	162
% 	160	439	162
% 	160	439	163
% 	161	438	162
% 	161	439	162
% 	161	439	163
% 	162	436	160
% 	162	436	161
% 	162	437	161
% 	163	436	160
% 	163	436	161
% 	163	437	161
% 	162	438	162
% 	162	439	162
% 	162	439	163
% 	163	438	162
% 	163	439	162
% 	163	439	163
% 	164	435	160
% 	165	435	160
% 	164	436	160
% 	164	436	161
% 	164	437	161
% 	165	436	160
% 	165	437	160
% 	165	437	161
% 	164	438	161
% 	165	438	161
% 	164	438	162
% 	164	439	162
% 	164	439	163
% 	165	438	162
% 	165	439	162
% 	166	436	160
% 	166	437	160
% 	166	437	161
% 	167	437	160
% 	166	438	161
% 	167	438	161
% 	167	439	161
% 	166	439	162
% 	167	439	162];
% 
% hold on gcf;
% scatter3(orig(:, 1), orig(:, 2), orig(:, 3), 30, 'filled', 'g');


