function check_abundances_sum(a_3D, plots_title, used_colormap)

% a_3D is a datacube with the abundances of each endmember; plots_title is
% the title of each plot; used_colormap is the employed colormap

[lines, columns, n_endmembers] = size(a_3D);

a_3D_sum = zeros(lines, columns);

for i = 1:lines
    for j = 1:columns
        for k = 1:n_endmembers
            
            a_3D_sum(i,j) = a_3D_sum(i,j) + a_3D(i,j,k);
            
        end
    end
end

plot_cube_3D(a_3D_sum, plots_title, used_colormap);