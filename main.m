clear all;
close all;
% Set the ICO level for high-resolution brain surface analysis. Default is 'ico4', 
% which uses a refined triangular mesh with 5,120 faces suitable for neuroimaging.
% This is a option in SBCI Pipeline.
icoLevel = 'ico4';

% Add required directories to the path for the initial setup.
addpath('./io');
addpath('./plotting');
addpath('./sfc');
addpath('./example_data/fsaverage_label/');
addpath('./example_data/SBCI_Individual_Subject_Outcome/');
addpath('./analysis');

% Load average data and surface from the SBCI pipeline.
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label', icoLevel);
sbci_surf = load_sbci_surface('example_data/fsaverage_label');

% Load functional and structural connectivity matrices.
smoothedFileName = sprintf('smoothed_sc_avg_0.005_%s.mat', icoLevel);
fcFileName = sprintf('fc_avg_%s.mat', icoLevel);
load(smoothedFileName);
load(fcFileName);

% Set and Display the atlas currently in use.
% atlas name can be found by `sbci_parc.atlas`
atlas_index = 42; % Example using 'aparc' atlas
disp(fprintf('the atlas current use is: %s', sbci_parc(atlas_index).atlas{1}));

% Define indices for regions to exclude from analysis (non-meaningful brain regions).
% Specific ROI Names can be found by `sbci_parc.names`
roi_exclusion_index = [1,36]; % Index 1: 'LH_missing', Index 36: 'RH_missing'


% Symmetrize the connectivity matrices & removing diagonal elements.
sc = sc + sc' - 2*diag(diag(sc));


% Transition to manipulation and visualization of connectivity matrices for brain network analysis.

% Visualize high-resolution FC and SC matrices.
plot_sbci_mat(fc, sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 1, 'clim', [-0.1, 0.35]);
title(['Continuous FC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 2, 'clim', [0, 3.5]);
title(['Continuous SC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

% Convert continuous SC and FC to ROI-based discrete matrices.
adjust_connectivity_res_for_fc(fc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);
adjust_connectivity_res_for_sc(sc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);

% Compute and display surface functional connectivity (SFC).
sfc_gbl = calculate_sfc_gbl(sc, fc, 'triangular', false);
sfc_loc = calculate_sfc_loc(sc, fc, sbci_parc(atlas_index), 'triangular', false);
plot_cortical_sfc(sfc_gbl, sfc_loc, sbci_surf, sbci_mapping, sbci_parc, atlas_index);

% Display any value on txt file to the surface.
plot_value_cortically(sbci_surf, sbci_mapping, 'example.txt')





