%Accumulate all the control points (not just the unique ones) into one cell
%array

nbr_levs = 8;

all_ctrlpts = cell(nbr_levs, 1);
for lev = 1:nbr_levs
    for cnt = 1:length(ctrl_pts_pointers{lev})
        all_ctrlpts{lev}(cnt) = reconstructed_control_points{lev}(ctrl_pts_pointers{lev}(cnt));
    end
    all_ctrlpts{lev} = all_ctrlpts{lev}';
end