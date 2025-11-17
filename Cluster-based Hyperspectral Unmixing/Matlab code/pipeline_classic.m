%% Application of the entire traditional pipeline on AVIRIS Cuprite dataset


clear all;
close all;
clc;

addpath(genpath(pwd));



%% Settings

data_mat_file = 'Cuprite97_data.mat';       % filename of the spectral datacube

pfa = 1e-5;                                 % False-alarm probability in NWHFC method

deterministic = 1;                          % 1 fixes the rng seed; 0 is for random seed

DEBUG = 1;                                  % 1 is for debug mode, 0 is without debug checks

used_colormap = 'jet';                      % other possible are 'gray', 'bone'

save = 0;

%% Set seed and other settings

if (deterministic)
    rng(1);
end



%% Import datacube

load(data_mat_file)

x_3D     = reshape_into_3D(x, Lines, Columns, DEBUG);

x_3D_HIP = hypercube(x_3D,wavlen(BANDS));             % Create hypercube object to perform useful operations
                                                      % using Hyperspectral Image Processing library


%% Plot RGB version of the cube

figure
imagesc(colorize(x_3D_HIP, 'Method', 'rgb', 'ContrastStretching', true))
axis equal
axis([1, Columns, 1, Lines])
title('Estimated RGB image of the area')


%% number of endmembers estimate

% HySime (The only used in this thesis analysis)

[w, Rn] = estNoise(x,'additive','off'); % Noise estimate
[kf, ~] = hysime(x, w, Rn, 'off');      % Estimated number of endmembers

disp(['The estimated N_endmembers through HySime is ', num2str(kf)])

% NWHFC (extra)

kf_NWHFC = countEndmembersHFC(x_3D_HIP, 'PFA', pfa, 'NoiseWhiten', 1);    % NWHFC is provided by Hyperspectral Image Processing library

disp(['The estimated N_endmembers through NWHFC is ', num2str(kf_NWHFC)])



%% Endmembers extraction and visual representation

% VCA algorithm

[Extract_endm, Select_pixels, x_proj] = VCA(x, 'Endmembers', kf, 'verbose', 'yes');

% Plot VCA-extracted endmembers as library

figure
colormap(used_colormap)

imagesc(Extract_endm)
colorbar
title('Endmembers signatures extracted through VCA (classic)')
xlabel('Endmember')
ylabel('Band')
xticks(1:kf)

% Plot VCA-extracted endmembers separately
N_endmembers = kf;
Extract_endm_extended          = NaN(size(wavlen,1), kf);
Extract_endm_extended(BANDS,:) = Extract_endm;

Elems_per_row = ceil(sqrt(N_endmembers));

figure

for (i = 1:N_endmembers)
    subplot(Elems_per_row,Elems_per_row, i)
    plot(wavlen, Extract_endm_extended(:,i), 'LineWidth',1)
    title(['(classic) Endmember ', num2str(i)])
    xlabel('Wavelength (nm)')
    ylabel('Reflectance (%)')
    xlim([min(wavlen), max(wavlen)])
    ylim([0 1])
    grid on
end


%% Original image and the selected pure pixels

% Retrieve the map of the pure pixels

sel_pixels_row = zeros(1,size(x,2));
sel_pixels_row(Select_pixels) = 1:N_endmembers;

sel_pixels_map = reshape_into_3D(sel_pixels_row, Lines, Columns, DEBUG);

% Plot the image

figure

imagesc(ones(250, 190))
colormap([1, 0.94, 0.86])
axis equal
axis([1, Columns, 1, Lines])
title('Purest pixels (classic)')
hold on

% Superimpose the positions of the purest pixels

for i = 1:Lines
    for j = 1:Columns
        if(sel_pixels_map(i,j) ~= 0)
            
            plot(j,i,'r.', 'MarkerSize', 15, 'Color', 'black')
            text(j+2,i,num2str(sel_pixels_map(i,j)), 'Color', 'black', 'FontSize',10)
            
            disp(['(',num2str(i),',',num2str(j),',',num2str(sel_pixels_map(i,j)),')'])
            
        end
    end
end


%% %% FCLS abundance estimation and visual representation

    
    % Abundance estimation and 3D-reshape
    a_FCLS = sunsal(Extract_endm, x, 'lambda',0, 'ADDONE','yes', 'POSITIVITY','yes', 'AL_iters',5000, 'TOL',1e-7, 'verbose','yes');
    
    a_FCLS_3D = reshape_into_3D(a_FCLS, Lines, Columns, DEBUG);
    
    % Plot abundances
    plot_cube_3D(a_FCLS_3D, 'FCLS (classic) - Abundance', used_colormap, [0 1]);
    
    
    FCLS_rec_error = abs(x - Extract_endm * a_FCLS);
    

    figure
    colormap(used_colormap)
    imagesc(FCLS_rec_error)
    colorbar
    title('FCLS - Reconstruction error per band per pixel (classic)')
    clim([0 0.16])
    xlabel('Pixel')
    ylabel('Band')
    
    FCLS_fro_error = norm(FCLS_rec_error, 'fro') / sqrt(numel(FCLS_rec_error));
    
    % Check violations of ANC
    check_ANC_violations(a_FCLS_3D, 'FCLS - ANC violations (classic) - Endm ', used_colormap);
    
    % Check violations of ASC
    check_ASC_violations(a_FCLS_3D, 'FCLS - ASC violations (classic)', used_colormap);
    title('FCLS - ASC violations (classic)')
    
    check_abundances_sum(a_FCLS_3D, 'FCLS - Abundances sum (classic)', used_colormap);
    title('FCLS - Abundances sum (classic)')

    
    % Compute total abundances
    N_pixels= Lines * Columns;
    all_abundances = plot_total_abundances(a_FCLS_3D, 'FCLS - Total endmember abundances (classic)');
    hold on
    yline(N_pixels/N_endmembers)
    hold off
  

    
    % Compute total number of active endmembers
    number_active_endmembers(a_FCLS_3D, 'FCLS - Total active endmembers (classic)', used_colormap);
    clim([0 16])
    title('FCLS - Total active endmembers (classic)')

    
    %% CLS abundance estimation


    % Abundance estimation and 3D-reshape
    a_CLS    = sunsal(Extract_endm, x, 'lambda',0, 'ADDONE','no', 'POSITIVITY','yes', 'AL_iters',2000, 'TOL',1e-7, 'verbose','yes');
    
    
    a_CLS_3D = reshape_into_3D(a_CLS, Lines, Columns, DEBUG);

    % Plot abundances
    plot_cube_3D(a_CLS_3D, 'CLS (classic) - Abundance', used_colormap, [0 1]);

    % Check reconstruction error

    CLS_rec_error = abs(x - Extract_endm * a_CLS);
   

    figure
    colormap(used_colormap)
    imagesc(CLS_rec_error)
    colorbar
    title('CLS - Reconstruction error per band per pixel (classic)')
    clim([0 0.16])
    xlabel('Pixel')
    ylabel('Band')
    
    CLS_fro_error = norm(CLS_rec_error, 'fro') / sqrt(numel(CLS_rec_error));

    %Check violations of ANC
   check_ANC_violations(a_CLS_3D, 'CLS - ANC violations (classic) - Endm', used_colormap);
   title('CLS - ANC violations (classic)')

    % Check violations of ASC
    check_ASC_violations(a_CLS_3D, 'CLS - ASC violations (classic)', used_colormap);
    title('CLS - ASC violations (classic)')
    
    check_abundances_sum(a_CLS_3D, 'CLS - Abundances sum (classic)', used_colormap);
    title('CLS - Abundances sum (classic)')


    % Compute total abundances

   all_CLS_abundances = plot_total_abundances(a_CLS_3D, 'CLS - Total endmember abundances (classic)');
   hold on
   yline(N_pixels/N_endmembers)
   hold off

    % Compute total number of active endmembers
    number_active_endmembers(a_CLS_3D, 'CLS - Total active endmembers (classic)', used_colormap);
    clim([0 16])
    title('CLS - Total active endmembers (classic)')

%% comparison boxplots

   figure
   hold on

   subplot(1,2,1)
   bar(sum(a_FCLS'), 'red')
   title('FCLS tot abundances (classic)')
  
   subplot(1,2,2)
   bar(sum(a_CLS'), 'red')
   title('CLS tot abundances (classic)')
   
   linkaxes()

%% preparing for matching


library_mat_file = 'USGS_spectra_LIB2_crop.mat';% filename of the spectral library
load(library_mat_file)
n_sel_spectra    = 6;    
n_endmembers = size(Extract_endm,2);
lib_card     = size(Library,2);



%% Compute the similarities of the library spectra with the extracted endmembers

% Initialize the matrix

simil_matrix = zeros(n_endmembers, lib_card);

% Perform the comparison spectrum-wise

for i = 1:n_endmembers
    if(DEBUG)
        disp(['Computing dissimilarity values for endmember ', num2str(i)]);
    end
    
    for j = 1:lib_card
        simil_matrix(i,j) = sam(Extract_endm(:,i), Library(BANDS,j));
    end
end

% Plot the similarities

figure
colormap(flipud(winter))
imagesc(simil_matrix)
colorbar
title('Dissimilarity values between the endmembers and the library spectra (classic)')
xlabel('Library spectrum')
ylabel('Endmember spectrum')
axis([1 lib_card 1 n_endmembers])



%% Select the most suitable library spectra for each extracted endmember

% Initialize the matrix with the indices of most similar library spectra

lib_sel_spectra = zeros(n_endmembers, lib_card);
lib_sel_angles  = zeros(n_endmembers, lib_card);

% Retrieve the most similar library spectra for each endmember

for (i = 1:n_endmembers)
    
    [sort_val, sort_ind] = sort(simil_matrix(i,:), 'ascend');
    
    lib_sel_spectra(i,:) = sort_ind;
    lib_sel_angles (i,:) = sort_val;
    
end

% Plot the endmembers and the identified library spectra

for (i = 1:n_endmembers)
    
    % Introductory operations for the plot
    
    figure
    hold on
    grid on
    
    legend_labels = {strcat(['Endmember ', num2str(i)])};

    % Plot the current spectrum

    curr_spectrum = ones(numel(wavlen),1);
    curr_spectrum(1:end) = NaN;
    curr_spectrum(BANDS) = Extract_endm(:,i);

    plot(wavlen, curr_spectrum, '-','LineWidth',1.75);
    
    % Plot the identified library spectra
    
    for (j = 1:n_sel_spectra)
        % Retrieve a suitable name in the legend for the library spectra
        
        curr_name = Lib_Names{lib_sel_spectra(i,j)};
        
        curr_inds = strfind(curr_name,'_');
        legend_labels{j+1} = [strrep(curr_name((curr_inds(2)+1):(curr_inds(end)-1)),'_',' '), ...
                              ' (', num2str(lib_sel_angles(i,j)), ')'];
        
        % Plot the spectra
        
        plot(wavlen, Library(:, lib_sel_spectra(i,j)), '-','LineWidth',0.75);
    end
    
    legend(legend_labels, 'Location', 'best')
    xlabel('Wavelength')
    ylabel('Reflectance')
    title(strcat(['(classic) Endmember ', num2str(i)]))
    ylim([0 1])
    xlim([min(wavlen), max(wavlen)])

end 

% Some useful prints

disp('The selected spectra are: ')
lib_sel_spectra(:,1:(n_sel_spectra*2))

disp('The respective angles are: ')
lib_sel_angles (:,1:(n_sel_spectra*2))


rmsSAE = norm(lib_sel_angles(:, 1), 'fro')/sqrt(kf)

%% Results 

disp(['FCLS frob error (without clusters) = ', num2str(FCLS_fro_error)]); % 0.005242
disp(['CLS frob error (without clusters) = ', num2str(CLS_fro_error)]); % 0.0046069
disp(['rmsSAE (without clusters) = ', num2str(rmsSAE)]); % 0.05996

%% Save all the plots

if save

    save_folder = "C:\Users\miche\Desktop\Michele\PoliMi\Magistrale\TESI\Codici personali\Immagini\Classic pipeline matlab";
    figHandles = findobj('Type', 'figure');

    for i = 1:length(figHandles)
        fig = figHandles(i);
        fig.WindowState = 'maximized';
    
        filename = fullfile(save_folder, sprintf('plot_%d.png', length(figHandles) - i + 1));
        exportgraphics(fig, filename, 'Resolution', 300);
        
        close(fig);
    end
end