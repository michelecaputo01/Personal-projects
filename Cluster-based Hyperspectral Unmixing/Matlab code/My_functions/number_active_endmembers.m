function number_active_endmembers(a_3D, plots_title, used_colormap)

% Plot the number of endmembers with nonzero abundances in each method

[lines columns n_endmembers] = size(a_3D);

active_endmembers       = (a_3D ~= 0);
total_active_endmembers = zeros(lines,columns);

for i = 1:lines
    for j = 1:columns
        for k = 1:n_endmembers
            
            total_active_endmembers(i,j) = total_active_endmembers(i,j) + active_endmembers(i,j,k);
            
        end
    end
end

plot_cube_3D(total_active_endmembers, plots_title, used_colormap); 