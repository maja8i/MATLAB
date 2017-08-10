%Accumulate all the control points (not just the unique ones) into one cell
%array
function all_ctrlpts = get_all_ctrlpts(control_points_array, ctrl_pts_pointers)

start_lvl = 3;
nbr_levs = 8;

all_ctrlpts = cell(nbr_levs, 1);
for lev = start_lvl:nbr_levs
    for cnt = 1:length(ctrl_pts_pointers{lev})
        all_ctrlpts{lev}(cnt) = control_points_array{lev}(ctrl_pts_pointers{lev}(cnt));
    end
    all_ctrlpts{lev} = all_ctrlpts{lev}';
end