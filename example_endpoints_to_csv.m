%% example_endpoints_to_csv.m
% Example script showing how to use endpoints_to_csv function
%
% MATLAB R2018a

%% Setup
clear all;
close all;
clc;

%% Add Required Paths
addpath(genpath(pwd));

%% Load SBCI Data
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label','ico4');

%% Set Parameters
% Choose grid file (determines coordinate system)
gridFile = './example_data/fsaverage_label/template_sphere_grid_ico4.mat';
% gridFile = './example_data/fsaverage_label/template_inflated_grid_ico4.mat';
% gridFile = './example_data/fsaverage_label/template_white_grid_ico4.mat';

% Input file
subjectFile = './example_data/SBCI_Individual_Subject_Outcome/mesh_intersections_ico4.mat';

% Select parcellation (e.g., Glasser = 25, CoCoNest-250 = 9)
parc_idx = 44;  % Change this to select different parcellation

% Output file
outputFile = './example_endpoints.csv';

%% Convert Endpoints to CSV File
endpoints_to_csv(gridFile, subjectFile, sbci_parc, sbci_mapping, parc_idx, outputFile, 'show_plot', true);

fprintf('\nDone! Output written to: %s\n', outputFile);
fprintf('CSV file contains columns:\n');
fprintf('  coord_in_x, coord_in_y, coord_in_z, surf_in, coord_out_x, coord_out_y, coord_out_z, surf_out, parc_in, parc_out\n');
fprintf('  (surf_in/surf_out: 0=left hemisphere, 1=right hemisphere)\n');

