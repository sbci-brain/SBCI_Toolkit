% MATLAB R2018a
% 
% get_endponit.m
% Script to calculate the cartesian coordinates of endpoints of streamlines connecting cortical surfaces
% at a given resolution and plot a subset of them onto the mesh. 

clear all;
close all;

addpath(genpath('./matlab'))

subject = 'sub-NDARINV0D5J9T8P/ses-2YearFollowUpYArm1/psc_sbci_final_files/sbci_connectome/';

% Load the grid of the given resolution: here we can choose from the
% following surfaces:
% 1. white surface - template_white_grid_ico4.mat
% 2. inflated white surface - template_inflated_grid_ico4.mat
% 3. sphere - template_sphere_grid_ico4.mat

ico_grid = load('./example_data/fsaverage_label/template_white_grid_ico4.mat');
%ico_grid = load('./endpoints_abcd_fixed/SBCI_AVG/template_sphere_grid_ico4.mat');

% Get the data for a single subject
subject_data = load(sprintf('./example_data/SBCI_Individual_Subject_Outcome/mesh_intersections_ico4.mat', subject));

% Get the cartesian coordinates for the given mesh
[coords_in, coords_out, surf_in, surf_out] = get_endpoint_coords(ico_grid, subject_data);

% Get a subset of endpoints, uniformly sampled
idx = datasample(1:size(coords_in,1), 100000); 

% Plot endpoints on left hemisphere 
left_ind = find(surf_in==0);
left_out = find(surf_out==0);
figure(1);clf;
trisurf(ico_grid.lh_T, ico_grid.lh_V(:,1), ico_grid.lh_V(:,2), ico_grid.lh_V(:,3), ...
    'FaceColor','flat', 'FaceVertexCData', ones(size(ico_grid.lh_V,1),3)*0.7, 'EdgeColor', 'none');
hold on;
scatter3(coords_in(left_ind,1), coords_in(left_ind,2),coords_in(left_ind,3),'.')
hold on;
scatter3(coords_out(left_out,1), coords_out(left_out,2),coords_out(left_out,3),'.')
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; 
material dull; shading interp;
view([1,0,0])


% % Plot the template mesh and endpoints in all brain
right_in = find(surf_in==1);
right_out = find(surf_out==1);
figure(2);clf;
trisurf(ico_grid.lh_T, ico_grid.rh_V(:,1), ico_grid.rh_V(:,2), ico_grid.rh_V(:,3), ...
    'FaceColor','flat', 'FaceVertexCData', ones(size(ico_grid.rh_V,1),3)*0.7, 'EdgeColor', 'none');
hold on;
scatter3(coords_in(right_in,1), coords_in(right_in,2),coords_in(right_in,3),'.')
hold on;
scatter3(coords_out(right_out,1), coords_out(right_out,2),coords_out(right_out,3),'.')
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; 
material dull; shading interp;
view([1,0,0])

% 
% % Make the plot look better
% daspect([1 1 1]); axis tight; camlight; axis vis3d off;
% lighting gouraud; 
% material dull; shading interp;
% view([1,0,0])
% 
% tiledlayout(4,1);
% ax1 = nexttile;
% plot(ax1,1:1149,cortical_ts(100,:));
% ax2 = nexttile;
% plot(ax2,1:1149,cortical_ts(101,:));
% ax3 = nexttile;
% plot(ax3,1:1149,cortical_ts(102,:));
% ax4 = nexttile;
% plot(ax4,1:1149,cortical_ts(103,:));

