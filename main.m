clear all;
close all;
% set ico level, 'ico4' as default
icoLevel = 'ico4';

% add path to local environment before first use
addpath('./io');
addpath('./plotting');
addpath('./sfc');
addpath('./example_data/');
addpath('./analysis');
% load sbci data from SBCI Pipeline result
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label', icoLevel);
sbci_surf = load_sbci_surface('example_data/fsaverage_label');
smoothedFileName = sprintf('smoothed_sc_avg_0.005_%s.mat', icoLevel);
fcFileName = sprintf('fc_avg_%s.mat', icoLevel);
load(smoothedFileName);
load(fcFileName);

% set specific atlas number that going to visualize
% atlas name can be found by `sbci_parc.atlas`
% display names for all parcellations
disp(strjoin(cellfun(@(c) c{1}, {sbci_parc.atlas}, 'UniformOutput', false), ', '));

atlas_index = 19; % aparc.a2005s in example

disp(fprintf('the atlas current use is: %s', sbci_parc(atlas_index).atlas{1}));




% Specify regions to be removed from the plot in the analysis
roi_mask_num = [1,36];

    
% 1. Display of High-Resolution SC/FC Data
%plot_connectivity(fc, sc, sbci_parc, atlas_num, roi_mask_num);

% 2. Convert continuous SC/FC to ROI based Discrete SC/FC
% adjust_connectivity_res(fc, sc, sbci_parc, atlas_index, sbci_mapping, roi_mask_num);

% 3. Compute SFC and Display it on the Cortical Surface
% calculate global and local SFC
sfc_gbl = calculate_sfc_gbl(sc, fc, 'triangular', true);
sfc_loc = calculate_sfc_loc(sc, fc, sbci_parc(atlas_index), 'triangular', true);
% Display on the Cortical Surface
plot_cortical_sfc(sfc_gbl, sfc_loc, sbci_surf, sbci_mapping, sbci_parc, atlas_index);

% 4. Display any value on txt file to the surfac
%plot_value_cortically(sbci_surf, sbci_mapping, 'example.txt')





