clear all;
close all;
% Set the ICO level for high-resolution brain surface analysis. Default is 'ico4', 
% which uses a refined triangular mesh with 5,120 faces suitable for neuroimaging.
% This is a option in SBCI Pipeline.
icoLevel = 'ico4';

% Add paths to the local environment before first use.
addpath('./io');
addpath('./plotting');
addpath('./sfc');
addpath('./example_data/');
addpath('./analysis');

% Load average data and surface from the SBCI pipeline.
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label', icoLevel);
sbci_surf = load_sbci_surface('example_data/fsaverage_label');

% load sbci individual data. These data can be found in 
% ‘example_data/SBCI_Individual_Subject_Outcome’ 
% The following two matrix are not included in the `example_data` folder, 
% please use the link to download the example data: 
% https://www.dropbox.com/scl/fo/35qxqnvo3xaq77yrpctz2/ABb2_GA0-dV0-7S1HmmTV4I?rlkey=c4koute3tpbvgzp2zt8i96unh&dl=0
smoothedFileName = sprintf('smoothed_sc_avg_0.005_%s.mat', icoLevel);
fcFileName = sprintf('fc_avg_%s.mat', icoLevel);
load(smoothedFileName);
load(fcFileName);

% set specific atlas for data visualization; 
% atlas name can be found by `sbci_parc.atlas`
atlas_index = 19; % aparc.a2005s in example
disp(fprintf('the atlas current use is: %s', sbci_parc(atlas_index).atlas{1}));




% Specify regions to be removed from the connectivity matrix. 
% Note that some brain regions is not meaningful in certain data analysis, 
% e.g., the corpus callosum in SC data analysis. 
% Specific ROI Names can be found by `sbci_parc.names`
roi_exclusion_index = [1,36]; 
% Currently, ROI 1 represents 'LH_G_cingulate-Isthmus'
% ROI 36 represents 'LH_Lat_Fissure-ant_sgt-ramus_horizontal'

% Symmetrizes the fc and sc matrix and clears the elements on the diagonals
fc = fc + fc' - 2*diag(diag(fc));
sc = sc + sc' - 2*diag(diag(sc));
% Normalize sc matrix
sc = sc/sum(sum(sc));

% TODO: write a summary sentence

% 1. Display of High-Resolution fc and sc matrix
plot_sbci_mat(fc, sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 1, 'clim', [-0.1, 0.35]);
title(['Continuous FC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 2, 'clim', [0, 3.5]);
title(['Continuous SC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

% 2. Convert continuous SC/FC to ROI based Discrete SC/FC
adjust_connectivity_res_for_fc(fc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);
adjust_connectivity_res_for_sc(sc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);

% 3. Compute SFC and display it on the fsaverage cortical surface
sfc_gbl = calculate_sfc_gbl(sc, fc, 'triangular', false);
sfc_loc = calculate_sfc_loc(sc, fc, sbci_parc(atlas_index), 'triangular', false);
% Display SFC the Cortical Surface
plot_cortical_sfc(sfc_gbl, sfc_loc, sbci_surf, sbci_mapping, sbci_parc, atlas_index);

% 4. Display any value on txt file to the surfac
%plot_value_cortically(sbci_surf, sbci_mapping, 'example.txt')





