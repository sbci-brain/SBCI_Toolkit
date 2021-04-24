clear all;
close all;

addpath('./example_data');
addpath('./io');
addpath('./plotting')

% load required SBCI data for mapping and analysis;
% This gives mapping from current resolution to a higher or lower
% resolution
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('example_data', 0.94);

% example to plot data on the inflated cortical surface
sbci_surf = load_sbci_surface('example_data');

%load data

load('./sfc_data/SBCI_SFC_global_rest1.mat')

for i=1:900
    curr_sfc = squeeze(sbci_sfc_gbl_tensor(:,i));
    plot_cortical(sbci_surf.inflated, sbci_mapping, curr_sfc,'figid',1);
    pause(0.1);
end