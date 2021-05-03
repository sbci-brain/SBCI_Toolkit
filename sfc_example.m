clear all;
close all;

addpath('./example_data');
addpath('./io');
addpath('./plotting');
addpath('./sfc');

% load required SBCI data for mapping and analysis;
% This gives mapping from current resolution to a higher or lower resolution
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('example_data', 0.94);

% load the surfaces to plot to
sbci_surf = load_sbci_surface('example_data');

% load the data
load('./example_data/example_fc.mat')
load('./example_data/example_sc.mat')

% calculate the different SFC measures
sfc_gbl = calculate_sfc_gbl(sc, fc, 'triangular', true);
sfc_loc = calculate_sfc_loc(sc, fc, sbci_parc(3), 'triangular', true);

% fill NaNs with 0s
sfc_gbl(isnan(sfc_gbl)) = 0;
sfc_loc(isnan(sfc_loc)) = 0;

% display SFC on the inflated surface
plot_cortical(sbci_surf.inflated, sbci_mapping, sfc_gbl, 'figid', 1, ...
    'legend', 'SFC global', 'clim', [0,1])

plot_cortical(sbci_surf.inflated, sbci_mapping, sfc_loc, 'figid', 2, ...
    'legend', 'SFC local (Desikan)', 'clim', [0,1])
