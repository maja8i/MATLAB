%------------------------------ User Inputs ------------------------------%

%Below assumes that the input point cloud PLY file is in the current 
%directory
ptcloud_name = 'frame_0000';    %Don't include ".ply" extension    
shifts = [1 2 3 4];
recoloured = 1; %0 => no; 1 => yes
shift_dir = 0;  %0 => x; 1 => y; 2 =>z

%-------------------------------------------------------------------------% 

%"Undistorted" point cloud that pc_error will compare to each time
ptcloud_in = [ptcloud_name '.ply'];

%Initialize array to store PSNR Y values
psnr_y = zeros(length(shifts), 1);
psnr_cntr = 1;

%Compute Y PSNR for the unshifted point cloud (i.e., original against 
%itself)
ptcloud_out = ptcloud_in;
%[~, sys_result] = system(['cd C:\MATLAB_Git-master\MPEG_Objective_Metrics_Input_Document &pc_error.exe -a ' ptcloud_in ' -b ' ptcloud_out ' -c -r 1023']);
%Below assumes that pc_error.exe is found in the current Matlab working
%directory
[~, sys_result] = system(['cd ' pwd ' &pc_error.exe -a ' ptcloud_in ' -b ' ptcloud_out ' -c -r 1023']);
start_psnr_str = strfind(sys_result, 'c[0],PSNRF') + 21;
end_psnr_str = strfind(sys_result, 'c[1],PSNRF') - 5;
psnr_y(psnr_cntr) = str2double(sys_result(start_psnr_str:end_psnr_str));
psnr_cntr = psnr_cntr + 1;

%Compute Y PSNR for each shifted point cloud against the original
for s = shifts
    %"Distorted" point cloud that pc_error will compare to
    if recoloured == 0
        if shift_dir == 0
            ptcloud_out = [ptcloud_name '_Xshift' num2str(s) '.ply'];
        elseif shift_dir == 1
            ptcloud_out = [ptcloud_name '_Yshift' num2str(s) '.ply'];
        elseif shift_dir == 2
            ptcloud_out = [ptcloud_name '_Zshift' num2str(s) '.ply'];
        end
    else
        if shift_dir == 0
            ptcloud_out = [ptcloud_name '_Xshift' num2str(s) '_recoloured.ply'];
        elseif shift_dir == 1
            ptcloud_out = [ptcloud_name '_Yshift' num2str(s) '_recoloured.ply'];
        elseif shift_dir == 2
            ptcloud_out = [ptcloud_name '_Zshift' num2str(s) '_recoloured.ply'];
        end
    end
    %[~, sys_result] = system(['cd C:\MATLAB_Git-master\MPEG_Objective_Metrics_Input_Document &pc_error.exe -a ' ptcloud_in ' -b ' ptcloud_out ' -c -r 1023']);
    %Below assumes that pc_error.exe is found in the current Matlab working
    %directory
    [~, sys_result] = system(['cd ' pwd ' &pc_error.exe -a ' ptcloud_in ' -b ' ptcloud_out ' -c -r 1023']);
    start_psnr_str = strfind(sys_result, 'c[0],PSNRF') + 21;
    end_psnr_str = strfind(sys_result, 'c[1],PSNRF') - 5;
    psnr_y(psnr_cntr) = str2double(sys_result(start_psnr_str:end_psnr_str));
    psnr_cntr = psnr_cntr + 1;
end

%Plot Y PSNR vs shifts in shift_dir direction
figure;
plot(shifts, psnr_y(2:end), '-o');
if strcmp('frame_0000', ptcloud_name)
    name_label = 'frame\_0000';
elseif strcmp('longdress_vox10_1330', ptcloud_name)
    name_label = 'longdress\_vox10\_1330';
end
set(gca, 'XTick', shifts);
if recoloured == 0
    if shift_dir == 0
        title({'Y PSNR vs Shift Along x Direction', ['for ' name_label]});
    elseif shift_dir == 1
        title({'Y PSNR vs Shift Along y Direction', ['for ' name_label]});
    elseif shift_dir == 2
        title({'Y PSNR vs Shift Along z Direction', ['for ' name_label]});
    end
else
    if shift_dir == 0
        title({'Y PSNR vs Shift Along x Direction', ['for Recoloured ' name_label]});
    elseif shift_dir == 1
        title({'Y PSNR vs Shift Along y Direction', ['for Recoloured ' name_label]});
    elseif shift_dir == 2
        title({'Y PSNR vs Shift Along z Direction', ['for Recoloured ' name_label]});
    end
end
xlabel('Shift Amount');
ylabel('Y PSNR (dB)');
grid on;
%Save plot figure in the current directory
if recoloured == 0
    if shift_dir == 0
        savefig([ptcloud_name '_Y_PSNR_vs_XShift']);
    elseif shift_dir == 1
        savefig([ptcloud_name '_Y_PSNR_vs_YShift']);
    elseif shift_dir == 2
        savefig([ptcloud_name '_Y_PSNR_vs_ZShift']);
    end
else
    if shift_dir == 0
        savefig([ptcloud_name '_Recoloured_Y_PSNR_vs_XShift']);
    elseif shift_dir == 1
        savefig([ptcloud_name '_Recoloured_Y_PSNR_vs_YShift']);
    elseif shift_dir == 2
        savefig([ptcloud_name '_Recoloured_Y_PSNR_vs_ZShift']);
    end
end



