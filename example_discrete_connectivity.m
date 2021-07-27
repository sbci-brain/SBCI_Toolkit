% load example data
load('./example_data/example_fc.mat')
load('./example_data/example_sc.mat')

% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('./example_data/', 0.94);
 
% convert into full matrices (for plotting)
fc = fc + fc' - 2*diag(diag(fc));
sc = sc + sc' - 2*diag(diag(sc));

% generate discrete connectivity 
dct_fc = parcellate_fc(fc, sbci_parc(3), sbci_mapping, 'roi_mask', [1,36]);
dct_sc = parcellate_sc(sc, sbci_parc(3), sbci_mapping, 'roi_mask', [1,36]);

%% Plot Continuous Connectivity
plot_sbci_mat(fc, sbci_parc(3), 'roi_mask', [1,36], 'figid', 1, 'clim', [-0.1, 0.35]);
title('Continuous FC (Desikan)')

plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(3), 'roi_mask', [1,36], 'figid', 2, 'clim', [0, 3.5]);
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
