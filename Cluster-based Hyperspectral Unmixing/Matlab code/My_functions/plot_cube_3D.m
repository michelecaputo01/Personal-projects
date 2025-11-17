function plot_cube_3D(a_3D, plots_title, used_colormap, lim)

% Plot a datacube, with each image associated to an endmember. a_3D has
% shape (lines x columns x n_endmembers); title is the string added at the
% top of each plot; used_colormap is the colormap of all the plots

% Preparatory data

[lines, columns, n_endmembers] = size(a_3D);

elems_per_row = ceil(sqrt(n_endmembers));

% Plot the images

figure
colormap(used_colormap)

for (i = 1:n_endmembers)

    subplot(elems_per_row,elems_per_row,i)
    imagesc(a_3D(:,:,i))
    if exist('lim', 'var')
        clim(lim)
    end
    axis equal
    title([plots_title, ' ', num2str(i)])
    axis([1, columns, 1, lines])
    colorbar

end