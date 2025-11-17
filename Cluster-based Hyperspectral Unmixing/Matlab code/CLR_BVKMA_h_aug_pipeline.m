%% Load data and cluster

clear all;
close all;
clc;

addpath(genpath(pwd));

save = 0;

data_mat_file = 'Cuprite97_data.mat';    

library_mat_file = 'USGS_spectra_LIB2_crop.mat'; % filename of the spectral library

clusterization = load("C:\Users\miche\Desktop\Michele\PoliMi\Magistrale\TESI\Codici personali\R code\clusterization_for_matlab_CLR.txt");

deterministic = 1;

DEBUG = 1;

used_colormap = 'jet'; 

if (deterministic)
    rng(1);
end

load(data_mat_file);

load(library_mat_file);

load("haug_outputs.mat")

x_3D     = reshape_into_3D(x, Lines, Columns, DEBUG);

x_3D_HIP = hypercube(x_3D,wavlen(BANDS)); 

cluster_map = reshape_into_3D(clusterization', Lines, Columns, DEBUG);

%% Algo for each cluster (separated)

k = max(clusterization);

VCA_pixels = [];

Extract_endm_total = [];

N_endmembers_for_cluster = zeros(1,k);

f_total_matched = zeros(size(haug_endm,2), size(x,2));

c_total_matched = zeros(size(haug_endm,2), size(x,2));

a_FCLS = [];

a_CLS = [];

rmsSAE = 0;

augmentation = [47, 53, 64, 128, 157, 197, 261, 310, 339, 348, 373, 393, 413, 499, 549, 573, 581 ]; %selected spectra of the library from ground-truth

for cl = 1:k

    indexes = find(clusterization == cl) ;

    local_x = x(:,indexes);

    x_aug = horzcat(local_x, Library(BANDS, augmentation));

    % Noise and N_endmembers estimation

    [w, Rn] = estNoise(x_aug,'additive','off');
    [kf, ~] = hysime(x_aug, w, Rn, 'off');

    disp(['The estimated N_endmembers in cluster ', num2str(cl) ,' through HySime is ', num2str(kf)]);

    % VCA algorithm

    [Extract_endm, Select_pixels, x_proj] = VCA(x_aug, 'Endmembers', kf, 'verbose', 'yes');

    cost_matrix = zeros(size(Extract_endm,2), size(haug_endm,2));

    for i = 1:size(Extract_endm,2)
        for j = 1:size(haug_endm,2)
            cost_matrix(i,j) = sam(Extract_endm(:,i), haug_endm(:,j));
        end
    end

    match = matchpairs(cost_matrix, 1000);

    Extract_endm_total = [Extract_endm_total , Extract_endm];

    Select_pixels(Select_pixels > length(indexes)) = NaN ;

    ind = nan(length(Select_pixels),1);

    for i = 1:length(ind)
        if ~isnan(Select_pixels(i))
            ind(i) = indexes(Select_pixels(i));
        end
    end

    VCA_pixels = [VCA_pixels ; ind];

    N_endmembers_for_cluster(cl) = kf;

    % Plot VCA-extracted endmembers separately
    
    N_endmembers = kf;
    Extract_endm_extended = NaN(size(wavlen,1), kf);
    Extract_endm_extended(BANDS,:) = Extract_endm;
    Elems_per_row = ceil(sqrt(size(haug_endm,2)));

    figure

    for i = 1:size(haug_endm,2)
        m = match(match(:,2)== i,1);
        if isempty(m)
            continue
        end
        subplot(Elems_per_row,Elems_per_row, i)
        plot(wavlen, Extract_endm_extended(:,m), 'LineWidth',1)
        title(['Endmember ', num2str(i), ' in cluster ', num2str(cl), ' (h-aug)'])
        xlabel('Wavelength (nm)')
        ylabel('Reflectance (%)')
        xlim([min(wavlen), max(wavlen)])
        ylim([0 1])
        grid on
    end

    % FCLS abundance estimation

    matrix_FCLS = sunsal(Extract_endm, local_x, 'lambda',0, 'ADDONE','yes', 'POSITIVITY','yes', 'AL_iters',5000, 'TOL',1e-7, 'verbose','yes');
    
    f = zeros(kf, Lines*Columns);

    f(:,indexes) = matrix_FCLS;

    a_FCLS = [a_FCLS ; f];

    f_3D = reshape_into_3D(f, Lines, Columns, DEBUG);

    check_ANC_violations(f_3D, ['FCLS (cluster ', num2str(cl) ') - ANC violations (h-aug) - Endm'], used_colormap);

    check_ASC_violations(f_3D, ['FCLS (cluster ', num2str(cl) ') - ASC violations (h-aug)'], used_colormap);
    title(['FCLS (cluster ', num2str(cl) ') - ASC violations (h-aug)'])

    check_abundances_sum(f_3D, ['FCLS (cluster ', num2str(cl) ') - Abundances sum (h-aug)'], used_colormap);
    title(['FCLS (cluster ', num2str(cl) ') - Abundances sum (h-aug)'])

    f_matched = zeros(size(haug_a_FCLS));

    for i = 1:size(match,1)
        f_matched(match(i,2),:) = f(match(i,1),:);
    end

    f_matched_3D = reshape_into_3D(f_matched, Lines, Columns, DEBUG);

    N_pixels = length(indexes);
    plot_total_abundances(f_matched_3D, ['FCLS (cluster ', num2str(cl) ') - Total abundances after endmembers matching (h-aug)']);
    hold on
    yline(N_pixels/N_endmembers)
    hold off

    f_total_matched = f_total_matched + f_matched ;

    number_active_endmembers(f_3D, ['FCLS (cluster ', num2str(cl) ') - Total active endmembers (h-aug)'], used_colormap);
    clim([0 16])
    title(['FCLS (cluster ', num2str(cl) ') - Total active endmembers (h-aug)'])
    
    % CLS abundance estimation

    matrix_CLS = sunsal(Extract_endm, local_x, 'lambda',0, 'ADDONE','no', 'POSITIVITY','yes', 'AL_iters',2000, 'TOL',1e-7, 'verbose','yes');

    c = zeros(kf, Lines*Columns);

    c(:,indexes) = matrix_CLS;

    a_CLS = [a_CLS ; c];

    c_3D = reshape_into_3D(c, Lines, Columns, DEBUG);

    check_ANC_violations(c_3D, ['CLS (cluster ', num2str(cl) ') - ANC violations (h-aug) - Endm'], used_colormap);

    check_ASC_violations(c_3D, ['CLS (cluster ', num2str(cl) ') - ASC violations (h-aug)'], used_colormap);
    title(['FCLS (cluster ', num2str(cl) ') - ASC violations (h-aug)'])

    check_abundances_sum(c_3D, ['CLS (cluster ', num2str(cl) ') - Abundances sum (h-aug)'], used_colormap);
    title(['CLS (cluster ', num2str(cl) ') - Abundances sum (h-aug)'])

    c_matched = zeros(size(haug_a_CLS));

    for i = 1:size(match,1)
        c_matched(match(i,2),:) = c(match(i,1),:);
    end

    c_matched_3D = reshape_into_3D(c_matched, Lines, Columns, DEBUG);

    N_pixels = length(indexes);
    plot_total_abundances(c_matched_3D, ['CLS (cluster ', num2str(cl) ') - Total abundances after endmembers matching (h-aug)']);
    hold on
    yline(N_pixels/N_endmembers)
    hold off

    c_total_matched = c_total_matched + c_matched ;

    number_active_endmembers(c_3D, ['CLS (cluster ', num2str(cl) ') - Total active endmembers (h-aug)'], used_colormap);
    clim([0 16])
    title(['CLS (cluster ', num2str(cl) ') - Total active endmembers (h-aug)'])

    % Comparison FCLS vs CLS

    figure
    hold on

    subplot(1,2,1)
    bar(sum(f_matched'), 'red')
    title(['FCLS total abundances (h-aug) - Cluster ', num2str(cl)])
  
    subplot(1,2,2)
    bar(sum(c_matched'), 'red')
    title(['CLS total abundances (h-aug) - Cluster ', num2str(cl)])
   
    linkaxes()
    hold off

    % Preparing libraries

    n_sel_spectra    = 3;    
    n_endmembers = size(Extract_endm,2);
    lib_card     = size(Library,2);

    % Compute similarities between endmembers and spectra

    simil_matrix = zeros(n_endmembers, lib_card);

    for i = 1:n_endmembers
        if(DEBUG)
        disp(['Computing dissimilarity values for endmember ', num2str(i)]);
        end
    
        for j = 1:lib_card
        simil_matrix(i,j) = sam(Extract_endm(:,i), Library(BANDS,j));
        end
    end

    figure
    colormap(flipud(winter))
    imagesc(simil_matrix)
    colorbar
    title(['Dissimilarity values between the endmembers (h-aug) and the library spectra in cluster ', num2str(cl)])
    xlabel('Library spectrum')
    ylabel('Endmember spectrum')
    axis([1 lib_card 1 n_endmembers])

    % Select the most similar spectra

    lib_sel_spectra = zeros(n_endmembers, lib_card);
    lib_sel_angles  = zeros(n_endmembers, lib_card);

    for i = 1:n_endmembers
    
        [sort_val, sort_ind] = sort(simil_matrix(i,:), 'ascend');
    
        lib_sel_spectra(i,:) = sort_ind;
        lib_sel_angles (i,:) = sort_val;
    
    end

    % Plot the identified spectra 

    for i = 1:n_endmembers
        figure
        hold on
        grid on

        m = find(match(:,1) == i);
        real_index = match(m,2);

        legend_labels = {strcat(['Endmember ', num2str(real_index)])};

        if isempty(real_index)
            legend_labels = {strcat('Endmember (not matched)')};
        end

        curr_spectrum = ones(numel(wavlen),1);
        curr_spectrum(1:end) = NaN;
        curr_spectrum(BANDS) = Extract_endm(:,i);

        plot(wavlen, curr_spectrum, '-','LineWidth',1.75);
    
        for j = 1:n_sel_spectra
            
            curr_name = Lib_Names{lib_sel_spectra(i,j)};
        
            curr_inds = strfind(curr_name,'_');
            legend_labels{j+1} = [strrep(curr_name((curr_inds(2)+1):(curr_inds(end)-1)),'_',' '),' (', num2str(lib_sel_angles(i,j)), ')'];
            
            plot(wavlen, Library(:, lib_sel_spectra(i,j)), '-','LineWidth',0.75);
        end
        
        legend(legend_labels, 'Location', 'best')
        xlabel('Wavelength')
        ylabel('Reflectance')
        if isempty(real_index)
            title(strcat(['Cluster ', num2str(cl),' (h-aug) - Endmember (not matched)']))
        else
            title(strcat(['Cluster ', num2str(cl),' (h-aug) - Endmember ', num2str(real_index)]))
        end
        xlim([min(wavlen), max(wavlen)])
        ylim([0 1])

    
    end 

    weight = length(indexes)/(Lines*Columns);
    
    rmsSAE = rmsSAE + weight * norm(lib_sel_angles(:, 1), 'fro')/sqrt(kf);

end 

%% Useful plots

% Location of the purest pixels

sel_pixels_row = zeros(1,size(x,2));
index_purest = [];

for i = 1:k
    n = N_endmembers_for_cluster(i);
    index_purest = [index_purest , 1:n];
end

for i = 1:length(VCA_pixels)
    if ~isnan(VCA_pixels(i)) 
        sel_pixels_row(VCA_pixels(i)) = index_purest(i);
    end
end

sel_pixels_map = reshape_into_3D(sel_pixels_row, Lines, Columns, DEBUG);

figure
imagesc(cluster_map)
colormap([1, 0, 0;
          1, 1, 0;
          0, 1, 0;
          0.5, 0.8, 1])
axis equal
axis([1, Columns, 1, Lines])
title('Purest pixels for each cluster (h-aug)')
hold on

for i = 1:Lines
    for j = 1:Columns
        if(sel_pixels_map(i,j) ~= 0)
            plot(j,i,'r.', 'MarkerSize', 15, 'Color', 'black')
            text(j+2,i,num2str(sel_pixels_map(i,j)), 'Color', 'black', 'FontSize', 10)
        end
    end
end

% FCLS visual representation

a_FCLS_3D = reshape_into_3D(a_FCLS, Lines, Columns, DEBUG);
number_active_endmembers(a_FCLS_3D, 'FCLS (grouped clusters) - Total active endmembers (h-aug)', used_colormap);
clim([0 16])
title('FCLS (grouped clusters) - Total active endmembers (h-aug)')

FCLS_rec_error = abs(x - Extract_endm_total * a_FCLS);
FCLS_fro_error = norm(FCLS_rec_error, 'fro') / sqrt(numel(FCLS_rec_error));

figure
colormap(used_colormap)
imagesc(FCLS_rec_error)
colorbar
title('FCLS - Reconstruction error per band per pixel (clustered h-aug)')
clim([0 0.16])
xlabel('Pixel')
ylabel('Band')

% CLS visual representation

a_CLS_3D = reshape_into_3D(a_CLS, Lines, Columns, DEBUG);
number_active_endmembers(a_CLS_3D, 'CLS (grouped clusters) - Total active endmembers (h-aug)', used_colormap);
clim([0 16])
title('CLS (grouped clusters) - Total active endmembers (h-aug)')

CLS_rec_error = abs(x - Extract_endm_total * a_CLS);
CLS_fro_error = norm(CLS_rec_error, 'fro') / sqrt(numel(CLS_rec_error));

figure
colormap(used_colormap)
imagesc(CLS_rec_error)
colorbar
title('CLS - Reconstruction error per band per pixel (clustered h-aug)')
clim([0 0.16])
xlabel('Pixel')
ylabel('Band')

% Total abundances after endmembers matching in FCLS 

f_total_matched_3D = reshape_into_3D(f_total_matched, Lines, Columns, DEBUG);

N_pixels = Lines*Columns;
plot_total_abundances(f_total_matched_3D, 'FCLS - Total abundances after endmembers matching (h-aug)');
hold on
yline(N_pixels/size(haug_endm,2))
hold off

% Total abundances after endmembers matching in CLS 

c_total_matched_3D = reshape_into_3D(c_total_matched, Lines, Columns, DEBUG);

N_pixels = Lines*Columns;
plot_total_abundances(c_total_matched_3D, 'CLS - Total abundances after endmembers matching (h-aug)');
hold on
yline(N_pixels/size(haug_endm,2))
hold off

% Total endmembers abundances on the map

lim = [0 1];

plot_cube_3D(f_total_matched_3D, 'FCLS (Matched) (h-aug) - Abund.', used_colormap, lim);

plot_cube_3D(c_total_matched_3D, 'CLS (Matched) (h-aug) - Abund.', used_colormap, lim);

% Difference of active endmembers in FCLS (classic vs clustered)

active_FCLS_haug = reshape(sum(haug_a_FCLS ~= 0), Lines, Columns);
active_FCLS_cluster = reshape(sum(a_FCLS ~= 0), Lines, Columns);

figure
imagesc(active_FCLS_haug - active_FCLS_cluster)
colorbar
axis equal
clim([-8 12])
title('FCLS - Difference of active endmembers (H-augmented - Clustered)')

% Difference of active endmembers in CLS (classic vs clustered)

active_CLS_haug = reshape(sum(haug_a_CLS ~= 0), Lines, Columns);
active_CLS_cluster = reshape(sum(a_CLS ~= 0), Lines, Columns);

figure
imagesc(active_CLS_haug - active_CLS_cluster)
colorbar
axis equal
clim([-8 12])
title('CLS - Difference of active endmembers (H-augmented - Clustered)')

% Quantitative results

disp(['FCLS frob error (with clusters) = ', num2str(FCLS_fro_error)]); % 0.0057294
disp(['CLS frob error (with clusters) = ', num2str(CLS_fro_error)]); % 0.004359
disp(['rmsSAE (with clusters) = ', num2str(rmsSAE)]); % 0.043665


%% Save all the plots

if save

    save_folder = "C:\Users\miche\Desktop\Michele\PoliMi\Magistrale\TESI\Codici personali\Immagini\CLR + BVKMA + h_aug pipeline matlab";
    figHandles = findobj('Type', 'figure');

    for i = 1:length(figHandles)
        fig = figHandles(i);
        fig.WindowState = 'maximized';
    
        filename = fullfile(save_folder, sprintf('plot_%d.png', length(figHandles) - i + 1));
        exportgraphics(fig, filename, 'Resolution', 300);
    
        close(fig);
    end
end
