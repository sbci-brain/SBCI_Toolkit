%%% This file demonstrates how to visualize a parcellation with the SBCI_Toolkit  

% Add SBCI_Toolkit files to path 
addpath(genpath(pwd));

% Load SBCI Data
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label','ico4');
sbci_surf = load_sbci_surface('example_data/fsaverage_label');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting Examples %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Each parcellation is stored in sbci_parc(i), and its corresponding 
% labels are stored in sbci_parc(i).labels 

% Plot Glasser parcellation on inflated brain surface
plot_parce(sbci_surf.inflated, sbci_mapping, sbci_parc(25).labels)

% A parcellation can also be plotted on the white or sphere surface
% Plot CoCoNest-250 parcellation on white surface
plot_parce(sbci_surf.white, sbci_mapping, sbci_parc(9).labels)

% A parcellation can also be outlined on top of another parcellation using
% the overlay_parc argument
% Outline Glasser parcellation on top of CoCoNest-250 parcellation
plot_parce(sbci_surf.inflated, sbci_mapping, sbci_parc(9).labels, ...
    'overlay_parc',sbci_parc(25))
