%% plot_endpoint.m
% Script to calculate the Cartesian coordinates of endpoints of 
% streamlines connecting cortical surfaces at a given resolution and 
% plot a subset of them onto the mesh.
%
% This script loads the necessary grid and subject data, computes the 
% endpoint coordinates, and then creates two beautiful figures showing 
% the endpoints on the left and right hemispheres, respectively.
%
% MATLAB R2018a

%% Setup
clear all;    % Clear workspace variables
close all;    % Close open figures
clc;          % Clear the command window

%% Add Required Paths
% Add directories that contain necessary functions and data.
addpath('io');
addpath('plotting');
addpath('sfc');
addpath('example_data/fsaverage_label/');
addpath('example_data/SBCI_Individual_Subject_Outcome/');
addpath('analysis');
addpath('endpoints');

%% Load Grid and Subject Data
% Choose one of the following templates (uncomment the one you wish to use):
% 1. White surface
%gridFile = './example_data/fsaverage_label/template_white_grid_ico4.mat';
% 2. Inflated white surface
% gridFile = './example_data/fsaverage_label/template_inflated_grid_ico4.mat';
% 3. Sphere surface
gridFile = './example_data/fsaverage_label/template_sphere_grid_ico4.mat';

ico_grid = load(gridFile);

% Load subject-specific streamline intersection data.
subjectFile = './example_data/SBCI_Individual_Subject_Outcome/mesh_intersections_ico4.mat';
subject_data = load(subjectFile);

%% Compute Endpoint Coordinates
% get_endpoint_coords returns:
%   - coords_in:  [Nx3] matrix of streamline starting coordinates
%   - coords_out: [Nx3] matrix of streamline ending coordinates
%   - surf_in:    [Nx1] indicator (0=left, 1=right) for starting point hemisphere
%   - surf_out:   [Nx1] indicator (0=left, 1=right) for ending point hemisphere
[coords_in, coords_out, surf_in, surf_out] = get_endpoint_coords(ico_grid, subject_data);

%% Optional: Subsample Endpoints for Plotting
% (Uncomment the following lines if you want to limit the number of endpoints
% plotted, which is useful when there are too many points.)
% max_plot_points = 100000;
% if size(coords_in,1) > max_plot_points
%     plot_idx = datasample(1:size(coords_in,1), max_plot_points, 'Replace', false);
% else
%     plot_idx = 1:size(coords_in,1);
% end

%% Plot Endpoints on the Left Hemisphere
% Identify indices for endpoints that belong to the left hemisphere.
left_in_idx  = find(surf_in == 0);
left_out_idx = find(surf_out == 0);

figure(1); clf;
set(gcf, 'Color', 'w');  % Set figure background to white

% Plot the left hemisphere mesh using trisurf.
trisurf(ico_grid.lh_T, ico_grid.lh_V(:,1), ico_grid.lh_V(:,2), ico_grid.lh_V(:,3), ...
    'FaceColor', 'flat', 'FaceVertexCData', ones(size(ico_grid.lh_V,1),3)*0.7, ...
    'EdgeColor', 'none');
hold on;

% Plot streamline start endpoints (in blue).
scatter3(coords_in(left_in_idx,1), coords_in(left_in_idx,2), coords_in(left_in_idx,3), ...
    20, 'b', 'filled');

% Plot streamline end endpoints (in red).
scatter3(coords_out(left_out_idx,1), coords_out(left_out_idx,2), coords_out(left_out_idx,3), ...
    20, 'r', 'filled');

% Enhance the plot appearance.
daspect([1 1 1]);      % Equal axis scaling
axis tight;           % Fit axes tightly around the data
camlight('headlight');% Add a light following the camera
lighting gouraud;     % Use smooth lighting
material dull;        % Set material properties
shading interp;       % Smooth color interpolation
view([1, 0, 0]);      % Set the camera view
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Left Hemisphere: Streamline Endpoints', 'FontSize', 12, 'FontWeight', 'bold');
legend('Mesh','Streamline Start','Streamline End', 'Location','best');
grid on;

%% Plot Endpoints on the Right Hemisphere
% Identify indices for endpoints that belong to the right hemisphere.
right_in_idx  = find(surf_in == 1);
right_out_idx = find(surf_out == 1);

figure(2); clf;
set(gcf, 'Color', 'w');  % White background for the figure

% Plot the right hemisphere mesh using trisurf.
% Note: Assuming the right hemisphere vertices are stored in rh_V.
trisurf(ico_grid.lh_T, ico_grid.rh_V(:,1), ico_grid.rh_V(:,2), ico_grid.rh_V(:,3), ...
    'FaceColor', 'flat', 'FaceVertexCData', ones(size(ico_grid.rh_V,1),3)*0.7, ...
    'EdgeColor', 'none');
hold on;

% Plot streamline start endpoints (in blue).
scatter3(coords_in(right_in_idx,1), coords_in(right_in_idx,2), coords_in(right_in_idx,3), ...
    20, 'b', 'filled');

% Plot streamline end endpoints (in red).
scatter3(coords_out(right_out_idx,1), coords_out(right_out_idx,2), coords_out(right_out_idx,3), ...
    20, 'r', 'filled');

% Enhance the plot appearance.
daspect([1 1 1]);
axis tight;
camlight('headlight');
lighting gouraud;
material dull;
shading interp;
view([1, 0, 0]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Right Hemisphere: Streamline Endpoints', 'FontSize', 12, 'FontWeight', 'bold');
legend('Mesh','Streamline Start','Streamline End', 'Location','best');
grid on;
