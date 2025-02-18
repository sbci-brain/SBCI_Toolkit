%% Connectivity Estimation via Riemannian Diffusion/Matern Kernel Smoothing
% This script computes continuous structural connectivity for a subject.
% The connectivity is estimated by smoothing the discrete connectivity 
% (derived from streamline endpoints) using a Riemannian diffusion (or Matern)
% kernel based on the Laplace–Beltrami eigenfunctions of the cortical surface.
%
% References:
% To Be Included


%% Initialization
clear;  % Clear workspace
close all;

% Add required directories to the path for the initial setup.
addpath('./io');
addpath('./plotting');
addpath('./sfc');
addpath('./example_data/fsaverage_label/');
addpath('./example_data/SBCI_Individual_Subject_Outcome/');
addpath('./concon_estimate');


% Define the number of nodes for left and right hemispheres (mesh resolution)
nndsL = 2562;
nndsR = 2562;
nnds = nndsL + nndsR;
icoLevel = 'ico4';

%% Load Laplace–Beltrami Eigenfunctions for Each Hemisphere
% Precomputed eigenfunctions (U) and eigenvalues (Lambda) for the left and
% right hemispheres are loaded from .mat files.
load('EV_LBO_ds_ico4_L.mat');  % Loads variables: Lambda_L, U_L, etc.
load('EV_LBO_ds_ico4_R.mat');  % Loads variables: Lambda_R, U_R, etc.

%% Set Kernel Parameters and Compute Kernel Matrices
% We define the kernel bandwidth (kappa) based on candidate values derived
% from the spectrum of the Laplace–Beltrami operator. These candidates span
% a range determined by the smallest and a mid-range eigenvalue.
nk = 10;  % Number of candidate kappa values

% Compute lower and upper bounds for kappa:
% kappa_min and kappa_max are derived using a diffusion kernel formulation.
kappa_min = sqrt(-2 * log(0.9) / Lambda_L(end));
kappa_max = sqrt(-2 * log(0.001) / Lambda_L(200));

% Generate a set of candidate kappa values (logarithmically spaced)
kappa_set = exp(linspace(log(kappa_min), log(kappa_max), nk))';

% Select a specific kappa value from the candidate set (e.g., 4th value)
kappa = kappa_set(4);

% Compute the Riemannian diffusion kernel matrices for the left and right hemispheres.
% These functions (compute_diffusion_kernel_matrix) must implement the kernel 
% computation based on the eigen-decomposition of the Laplace–Beltrami operator.
KM_L = compute_diffusion_kernel_matrix(kappa, Lambda_L, U_L);
KM_R = compute_diffusion_kernel_matrix(kappa, Lambda_R, U_R);

% Alternatively, to use a Riemannian Matern kernel, uncomment the following lines:
% nu = 3;  % Matern smoothness parameter (possible values: 1, 2, 3)
% KM_L = compute_matern_kernel_matrix(nu, kappa, Lambda_L, U_L);
% KM_R = compute_matern_kernel_matrix(nu, kappa, Lambda_R, U_R);

%% Load Streamline Data and Construct Adjacency Matrices
% The file 'mesh_intersections_ico4.mat' contains streamline endpoint data:
%   - surf_in: indicator for hemisphere of the streamline's starting point (0 = left, 1 = right)
%   - surf_out: indicator for hemisphere of the streamline's ending point
%   - vtx_in: vertex index of the streamline's starting point
%   - vtx_out: vertex index of the streamline's ending point
load('mesh_intersections_ico4.mat', 'surf_in', 'surf_out', 'vtx_in', 'vtx_out');

% --- Left-to-Left Connectivity ---
% Identify streamlines that start and end in the left hemisphere.
ind = (surf_in == 0) & (surf_out == 0);
A11 = accumarray([vtx_in(ind)', vtx_out(ind)'], 1, [nndsL, nndsL]);
A11 = A11 + A11';  % Make the matrix symmetric

% --- Right-to-Right Connectivity ---
% Identify streamlines that start and end in the right hemisphere.
ind = (surf_in == 1) & (surf_out == 1);
A22 = accumarray([vtx_in(ind)', vtx_out(ind)'], 1, [nndsR, nndsR]);
A22 = A22 + A22';

% --- Inter-Hemispheric Connectivity (Left-to-Right) ---
% Identify streamlines connecting the left and right hemispheres.
ind = (surf_in == 0) & (surf_out == 1);
A12 = accumarray([vtx_in(ind)', vtx_out(ind)'], 1, [nndsL, nndsR]);

% Identify streamlines connecting right-to-left and add to A12 (transpose later)
ind = (surf_in == 1) & (surf_out == 0);
A21 = accumarray([vtx_in(ind)', vtx_out(ind)'], 1, [nndsR, nndsL]);
A12 = A12 + A21';  % Combine the two directions
clear A21;

%% Compute the Discretized Smoothed Connectivity (Density) Matrix
% The discrete connectivity (adjacency) matrices are smoothed by sandwiching 
% them between the kernel matrices. This yields intensity matrices that represent
% the continuous connectivity.

% --- Compute Intensity Matrices ---
% Left hemisphere:
IM11 = KM_L * A11 * KM_L;
% Right hemisphere:
IM22 = KM_R * A22 * KM_R;
% Inter-hemispheric connectivity:
IM12 = KM_L * A12 * KM_R;
% Note: IM21 is the transpose of IM12.

% --- Normalize Intensities to Obtain Densities ---
% Sum only the positive intensity values to compute the normalizing constant.
total_IM = sum(IM11(IM11 > 0)) + sum(IM22(IM22 > 0)) + (2 * sum(IM12(IM12 > 0)));

% Convert the intensity matrices to density matrices by normalizing.
DM11 = (IM11 ./ total_IM) .* (IM11 > 0);
DM12 = (IM12 ./ total_IM) .* (IM12 > 0);
DM22 = (IM22 ./ total_IM) .* (IM22 > 0);
% DM21 is simply DM12 transposed.

% Clear the intermediate intensity matrices.
clear IM11 IM12 IM22;

%% Assemble the Full Discretized Density Matrix
% Create a full connectivity (density) matrix DM that includes:
%   - Left-to-left connectivity in the top-left block,
%   - Right-to-right connectivity in the bottom-right block,
%   - Inter-hemispheric connectivity in the off-diagonal blocks.
DM = zeros(nnds, nnds);
DM(1:nndsL, 1:nndsL) = DM11;
DM(nndsL+1:end, nndsL+1:end) = DM22;
DM(1:nndsL, nndsL+1:end) = DM12;
DM(nndsL+1:end, 1:nndsL) = DM12';  % Ensure symmetry

% The output matrix DM is a discretized density matrix representing the 
% smoothed, continuous structural connectivity of the subject.


%plot DM
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label', icoLevel);
sbci_surf = load_sbci_surface('example_data/fsaverage_label');

% Set and Display the atlas currently in use.
% atlas name can be found by `sbci_parc.atlas`
atlas_index = 42; % Example using 'aparc' atlas
disp(fprintf('the atlas current use is: %s', sbci_parc(atlas_index).atlas{1}));

% Define indices for regions to exclude from analysis (non-meaningful brain regions).
% Specific ROI Names can be found by `sbci_parc.names`
roi_exclusion_index = [1,36]; % Index 1: 'LH_missing', Index 36: 'RH_missing'


% Visualize DM with node-reorder based on the sbci_parc(atlas_index).atlas{1}
plot_sbci_mat(log(10^7*DM+1), sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 1,'clim', [0, 3.5]);
title(['RDK Estimated ConCon (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

