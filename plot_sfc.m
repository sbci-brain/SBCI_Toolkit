clear all;
close all;

addpath('./example_data');
addpath('./io');
addpath('./plotting')

% load required SBCI data for mapping and analysis;
% This gives mapping from current resolution to a higher or lower
% resolution
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('example_data', 0.94);

sbci_surf = load_sbci_surface('example_data');


%load Structural and Functional Coupling and plot it
load example_sfc.mat
plot_cortical(sbci_surf.inflated, sbci_mapping, sfc);