%% parce_endpoint_overlay.m
% Script to overlay streamline endpoints on a parcellation plot with transparency.
% Allows filtering endpoints to show only those connecting selected parcel pairs.
%
% MATLAB R2018a

%% Setup
clear all;    % Clear workspace variables
close all;    % Close open figures
clc;          % Clear the command window

%% Add Required Paths
addpath(genpath(pwd));

%% Load SBCI Data
[sbci_parc, sbci_mapping, ~] = load_sbci_data('example_data/fsaverage_label','ico4');
sbci_surf = load_sbci_surface('example_data/fsaverage_label');

%% Load Grid and Subject Data
% Choose surface type (should match the surface used for parcellation)
%gridFile = './example_data/fsaverage_label/template_white_grid_ico4.mat';
gridFile = './example_data/fsaverage_label/template_inflated_grid_ico4.mat';

ico_grid = load(gridFile);

% Load subject-specific streamline intersection data
subjectFile = './example_data/SBCI_Individual_Subject_Outcome/mesh_intersections_ico4.mat';
subject_data = load(subjectFile);

%% Compute Endpoint Coordinates
[coords_in, coords_out, ~, ~] = get_endpoint_coords(ico_grid, subject_data);

%% Map Endpoints to Parcels
% Get node indices from intersection data
% Note: vtx_in and vtx_out are local hemisphere indices (1-2562 for each hemisphere)
% surf_in and surf_out indicate hemisphere (0=left, 1=right)
% parc_labels is a global array covering both hemispheres (1-5124 total)
vtx_in = subject_data.vtx_in;
vtx_out = subject_data.vtx_out;
surf_in = subject_data.surf_in;
surf_out = subject_data.surf_out;

% Select which parcellation to use (e.g., Glasser = 25, CoCoNest-250 = 9)
parc_idx = 44;  % Change this to select different parcellation
parc_labels = sbci_parc(parc_idx).labels;

% Get number of left hemisphere vertices
n_lh = double(sbci_mapping.shape(2));  % Number of left hemisphere vertices (typically 2562 for ico4)

% Convert local hemisphere indices to global indices
% Left hemisphere (surf_in == 0): global_idx = vtx_in (1 to n_lh)
% Right hemisphere (surf_in == 1): global_idx = n_lh + vtx_in (n_lh+1 to 2*n_lh)
% Convert to double to avoid integer type mismatch
global_idx_in = double(vtx_in);
global_idx_in(surf_in == 1) = n_lh + double(vtx_in(surf_in == 1));

global_idx_out = double(vtx_out);
global_idx_out(surf_out == 1) = n_lh + double(vtx_out(surf_out == 1));

% Check if global indices are within valid range
max_node_idx = length(parc_labels);
valid_in = (global_idx_in > 0) & (global_idx_in <= max_node_idx);
valid_out = (global_idx_out > 0) & (global_idx_out <= max_node_idx);

% Map global indices to parcel labels
parc_in = zeros(size(vtx_in));
parc_out = zeros(size(vtx_out));

parc_in(valid_in) = parc_labels(global_idx_in(valid_in));
parc_out(valid_out) = parc_labels(global_idx_out(valid_out));

% Debug: Print some sample mappings to verify (before filtering)
fprintf('Sample endpoint-to-parcel mapping (before filtering):\n');
if ~isempty(vtx_in)
    sample_idx = min(10, length(vtx_in));
    fprintf('  First %d endpoints:\n', sample_idx);
    for i = 1:sample_idx
        if valid_in(i) && valid_out(i)
            fprintf('    Endpoint %d: vtx_in=%d -> parcel %d, vtx_out=%d -> parcel %d\n', ...
                i, vtx_in(i), parc_in(i), vtx_out(i), parc_out(i));
        end
    end
end

%% Select Parcel Pairs to Display
% Option 1: Show all endpoints (set to empty)
%selected_pairs = [23,16; 10,16; 16,6;6,14];  % Empty means show all
selected_pairs = [12,33;19,15;28,15;47,68];  % Empty means show all
%selected_pairs = [1,2];  % Empty means show all
% Option 2: Show endpoints connecting specific parcel pairs
% Format: [parcel1, parcel2] - each row is a pair
% Example: Show connections between parcels 1 and 2, and parcels 3 and 4
% selected_pairs = [1, 2; 3, 4];

% Option 3: Show endpoints connecting one parcel to all others
% target_parcel = 1;
% all_parcels = unique(parc_labels);
% other_parcels = setdiff(all_parcels, target_parcel);
% selected_pairs = [repmat(target_parcel, length(other_parcels), 1), other_parcels(:)];

%% Extract Selected Parcels and Filter Endpoints
if isempty(selected_pairs)
    % Show all endpoints and all parcels
    selected_parcels = [];
    plot_idx = 1:length(vtx_in);
    fprintf('Showing all %d endpoints and all parcels\n', length(plot_idx));
else
    % Extract unique parcel IDs from selected pairs
    selected_parcels = unique(selected_pairs(:));
    fprintf('Selected parcels to display: %s\n', mat2str(selected_parcels));
    
    % Filter to show only selected parcel pairs
    plot_mask = false(size(vtx_in));
    for i = 1:size(selected_pairs, 1)
        pair = selected_pairs(i, :);
        % Check both directions (A->B and B->A)
        match = (parc_in == pair(1) & parc_out == pair(2)) | ...
                (parc_in == pair(2) & parc_out == pair(1));
        plot_mask = plot_mask | match;
    end
    plot_idx = find(plot_mask);
    
    fprintf('Filtered to %d endpoints connecting selected parcel pairs\n', length(plot_idx));
    if isempty(plot_idx)
        warning('No endpoints found for selected parcel pairs. Showing all endpoints instead.');
        plot_idx = 1:length(vtx_in);
    else
        % Debug: Show examples of selected streamlines and their parcel connections
        fprintf('\n=== Selected Streamlines Information ===\n');
        
        % Show 10 examples of selected streamlines
        filtered_parc_in = parc_in(plot_idx);
        filtered_parc_out = parc_out(plot_idx);
        sample_size = min(10, length(plot_idx));
        fprintf('Example %d selected streamlines:\n', sample_size);
        for i = 1:sample_size
            idx = plot_idx(i);
            fprintf('  Streamline %d: Parcel %d <-> Parcel %d', ...
                i, filtered_parc_in(i), filtered_parc_out(i));
            if valid_in(idx) && valid_out(idx)
                fprintf(' (vtx_in=%d, vtx_out=%d)', vtx_in(idx), vtx_out(idx));
            end
            fprintf('\n');
        end
        
        % Show count for each selected parcel pair
        fprintf('\nStreamline counts for each selected parcel pair:\n');
        for i = 1:size(selected_pairs, 1)
            pair = selected_pairs(i, :);
            % Count streamlines connecting this pair (both directions)
            count = sum((filtered_parc_in == pair(1) & filtered_parc_out == pair(2)) | ...
                        (filtered_parc_in == pair(2) & filtered_parc_out == pair(1)));
            fprintf('  Parcel pair [%d, %d]: %d streamlines\n', pair(1), pair(2), count);
        end
        
        % Also show all unique pairs found (in case there are more than selected)
        unique_pairs_found = unique([filtered_parc_in, filtered_parc_out], 'rows');
        if size(unique_pairs_found, 1) > size(selected_pairs, 1)
            fprintf('\nNote: Additional parcel pairs found in filtered endpoints:\n');
            for i = 1:size(unique_pairs_found, 1)
                pair = unique_pairs_found(i, :);
                % Check if this pair is in selected_pairs
                is_selected = any((selected_pairs(:,1) == pair(1) & selected_pairs(:,2) == pair(2)) | ...
                                 (selected_pairs(:,1) == pair(2) & selected_pairs(:,2) == pair(1)));
                if ~is_selected
                    count = sum((filtered_parc_in == pair(1) & filtered_parc_out == pair(2)) | ...
                                (filtered_parc_in == pair(2) & filtered_parc_out == pair(1)));
                    fprintf('  Parcel pair [%d, %d]: %d streamlines\n', pair(1), pair(2), count);
                end
            end
        end
        fprintf('==========================================\n\n');
    end
end

%% Optional: Subsample if too many points
max_plot_points = 50000;  % Maximum number of endpoints to plot
if length(plot_idx) > max_plot_points
    plot_idx = datasample(plot_idx, max_plot_points, 'Replace', false);
    fprintf('Subsampled to %d endpoints for visualization\n', length(plot_idx));
end

%% Plot Parcellation with Transparency and Overlay Endpoints
% Set transparency level (0 = fully transparent, 1 = opaque)
parcel_alpha = 0.1;  % Adjust this value to control transparency

% Select surface type (should match gridFile)
surface_type = 'inflated';  % Options: 'white', 'inflated', 'sphere'

% Don't mask parcels - show all parcels, but highlight selected ones
% Pass selected_parcels to plotting function so it can color them differently

% Plot with transparency and endpoints
% Note: coords_in, coords_out, surf_in, surf_out are already filtered by plot_idx
plot_parce_with_endpoints(sbci_surf.(surface_type), sbci_mapping, parc_labels, ...
    coords_in(plot_idx, :), coords_out(plot_idx, :), ...
    surf_in(plot_idx), surf_out(plot_idx), ...
    'parcel_alpha', parcel_alpha, ...
    'my_cols', [], ...  % Use default colors, or provide custom colormap
    'endpoint_size', 20, ...  % Size of endpoint markers
    'show_lines', false, ...  % Set to true to draw lines connecting endpoint pairs
    'selected_parcels', selected_parcels);  % Pass selected parcels to color them differently

fprintf('Plotting complete!\n');
fprintf('Tip: Adjust parcel_alpha to change transparency (current: %.2f)\n', parcel_alpha);
fprintf('Tip: Modify selected_pairs to filter endpoints by parcel connections\n');

