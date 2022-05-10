clear all;
close all;

addpath(genpath('../SBCI_Toolkit'))

% load example data
load('./example_data/example_fc.mat')
load('./example_data/example_sc.mat')

% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('./example_data/', 0.94);
 
% convert into full matrices (for plotting)
fc = fc + fc' - 2*diag(diag(fc));
sc = sc + sc' - 2*diag(diag(sc));

% In this example, we use Desikan atlas. The 'aparc' atlas in sbci_parc corresponds to the Desikan atlas and
% sbci_parc(11);

% We use 'parcellate_fc'/'parcellate_sc' to calculate discrete FC/SC. 
% Note that here the fourth argument in these functions is 'roi_mask',
% which is a vector of ROI labels to remove;

my_ROI_list = sbci_parc(11).names; %ROI names

%for Desikan, ROIs 1 and 36 are 'LH_missing' and 'RH_missing'. So let's remove them 
my_ROI_list([1,36],:) = []; %updated ROI names

dct_fc = parcellate_fc(fc, sbci_parc(11), sbci_mapping, 'roi_mask', [1,36]);
dct_sc = parcellate_sc(sc, sbci_parc(11), sbci_mapping, 'roi_mask', [1,36]);

%% Plot Continuous Connectivity
plot_sbci_mat(fc, sbci_parc(11), 'roi_mask', [1,36], 'figid', 1, 'clim', [-0.1, 0.35]);
title('Continuous FC (Desikan)')

plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(11), 'roi_mask', [1,36], 'figid', 2, 'clim', [0, 3.5]);
title('Continuous SC (Desikan)')

%% Plot Discrete FC
figure(3);
imagesc(dct_fc);

% remove ticks
xticks([]); xticklabels([]);
yticks([]); yticklabels([]);

% colours and aspect ratio
axis square;
daspect([1,1,1]);
set(gca, 'CLim', [-0.1, 0.35]);
colorbar();

title('Discrete FC (Desikan)')

%% Plot Discrete SC
figure(4);
imagesc(log((10^7*dct_sc) + 1));

% remove ticks
xticks([]); xticklabels([]);
yticks([]); yticklabels([]);

% colours and aspect ratio
axis square;
daspect([1,1,1]);
set(gca, 'CLim', [0, 3.5]);
colorbar();

title('Discrete SC (Desikan)')
