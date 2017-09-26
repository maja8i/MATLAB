%Accumulate the control points for ALL the occupied octree cell corners at
%each level (not just the control points for the unique corners) into one 
%cell array
function all_ctrlpts = get_all_ctrlpts(control_points_array, ctrl_pts_pointers, start_lvl, max_lvl)

all_ctrlpts = cell(max_lvl, 1);
for lev = start_lvl:max_lvl
%     for cnt = 1:length(ctrl_pts_pointers{lev})
%         all_ctrlpts{lev}(cnt) = control_points_array{lev}(ctrl_pts_pointers{lev}(cnt));
%     end
    all_ctrlpts{lev} = control_points_array{lev}(ctrl_pts_pointers{lev});
%     all_ctrlpts{lev} = all_ctrlpts{lev}';
end