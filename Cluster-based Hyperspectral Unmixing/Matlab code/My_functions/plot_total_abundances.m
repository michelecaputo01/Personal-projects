function total_abundances = plot_total_abundances(a_3D, plots_title)

% Plot a barplot with the total abundances of each endmember along the
% entire image; a_3D is a datacube with abundances for each endmember

% Preparatory data

n_endmembers     = size(a_3D, 3);
total_abundances = zeros(n_endmembers,1);

% Gather the abundances

for i = 1:n_endmembers
    total_abundances(i) = sum(sum(a_3D(:,:,i)));
end

figure
bar(1:n_endmembers, total_abundances, 'red');
title(plots_title)
xlabel('Endmember')
ylabel('Total abundance')
xticks(1:n_endmembers)