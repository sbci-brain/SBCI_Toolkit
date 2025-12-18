% MATLAB R2018a
%
% FUNCTION NAME:
%   endpoints_to_csv
%
% DESCRIPTION:
%   Convert streamline endpoints from mesh_intersections to CSV file format.
%   Each line contains: coord_in(3), surf_in, coord_out(3), surf_out, parc_in, parc_out
%
% INPUT:
%   gridFile - (string) Path to grid file (e.g., 'template_sphere_grid_ico4.mat')
%   subjectFile - (string) Path to mesh intersections file (e.g., 'mesh_intersections_ico4.mat')
%   sbci_parc - (struct array) Parcellation structures from load_sbci_data
%   sbci_mapping - (struct) Mapping structure from load_sbci_data
%   parc_idx - (scalar) Index of parcellation to use
%   outputFile - (string) Path to output CSV file
%   varargin - Optional arguments:
%       'show_plot' - (logical) Whether to display visualization plot (default: false)
%
% OUTPUT:
%   None (writes to file)
%
% ASSUMPTIONS AND LIMITATIONS:
%   The grid file, subject file, and parcellation must all be from the same resolution (ico4).

function endpoints_to_csv(gridFile, subjectFile, sbci_parc, sbci_mapping, parc_idx, outputFile, varargin)

% Parse optional arguments
p = inputParser;
addParameter(p, 'show_plot', false, @islogical);
parse(p, varargin{:});
params = p.Results;

%% Load Grid and Subject Data
ico_grid = load(gridFile);
subject_data = load(subjectFile);

%% Compute Endpoint Coordinates
fprintf('Computing endpoint coordinates...\n');
[coords_in, coords_out, ~, ~] = get_endpoint_coords(ico_grid, subject_data);

%% Get Vertex Indices and Hemisphere Information
vtx_in = subject_data.vtx_in;
vtx_out = subject_data.vtx_out;
surf_in = subject_data.surf_in;  % 0=left, 1=right
surf_out = subject_data.surf_out;  % 0=left, 1=right

%% Map Endpoints to Parcels
fprintf('Mapping endpoints to parcels...\n');
parc_labels = sbci_parc(parc_idx).labels;

% Get number of left hemisphere vertices
n_lh = double(sbci_mapping.shape(2));  % Number of left hemisphere vertices (typically 2562 for ico4)

% Convert local hemisphere indices to global indices
% Left hemisphere (surf_in == 0): global_idx = vtx_in (1 to n_lh)
% Right hemisphere (surf_in == 1): global_idx = n_lh + vtx_in (n_lh+1 to 2*n_lh)
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

%% Visualize Endpoints on Hemispheres
if params.show_plot
    fprintf('Creating visualization plot...\n');
    try
    % Subsample to 10% of endpoints for visualization
    num_total = length(vtx_in);
    num_plot = max(1, round(0.1 * num_total));
    plot_indices = datasample(1:num_total, num_plot, 'Replace', false);
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    
    % Left Hemisphere
    subplot(1, 2, 1);
    trisurf(ico_grid.lh_T, ico_grid.lh_V(:,1), ico_grid.lh_V(:,2), ico_grid.lh_V(:,3), ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    left_in = intersect(find(surf_in == 0), plot_indices);
    left_out = intersect(find(surf_out == 0), plot_indices);
    if ~isempty(left_in)
        scatter3(coords_in(left_in, 1), coords_in(left_in, 2), coords_in(left_in, 3), ...
            15, [0 0 0], 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 0.3);
    end
    if ~isempty(left_out)
        scatter3(coords_out(left_out, 1), coords_out(left_out, 2), coords_out(left_out, 3), ...
            15, [0 0 0], 'filled', 'MarkerEdgeColor', 'b', 'LineWidth', 0.3);
    end
    daspect([1 1 1]); axis off; axis tight; view(90, 0);
    title(sprintf('Left Hemisphere\n%d endpoints (showing %d)', ...
        sum(surf_in == 0) + sum(surf_out == 0), length(left_in) + length(left_out)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    camlight('headlight'); lighting gouraud;
    
    % Right Hemisphere
    subplot(1, 2, 2);
    trisurf(ico_grid.rh_T, ico_grid.rh_V(:,1), ico_grid.rh_V(:,2), ico_grid.rh_V(:,3), ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    right_in = intersect(find(surf_in == 1), plot_indices);
    right_out = intersect(find(surf_out == 1), plot_indices);
    if ~isempty(right_in)
        scatter3(coords_in(right_in, 1), coords_in(right_in, 2), coords_in(right_in, 3), ...
            15, [0 0 0], 'filled', 'MarkerEdgeColor', 'r', 'LineWidth', 0.3);
    end
    if ~isempty(right_out)
        scatter3(coords_out(right_out, 1), coords_out(right_out, 2), coords_out(right_out, 3), ...
            15, [0 0 0], 'filled', 'MarkerEdgeColor', 'b', 'LineWidth', 0.3);
    end
    daspect([1 1 1]); axis off; axis tight; view(-90, 0);
    title(sprintf('Right Hemisphere\n%d endpoints (showing %d)', ...
        sum(surf_in == 1) + sum(surf_out == 1), length(right_in) + length(right_out)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    camlight('headlight'); lighting gouraud;
    
    sgtitle(sprintf('Extracted Endpoints Visualization (showing 10%%: %d/%d)', num_plot, num_total), ...
        'FontSize', 14, 'FontWeight', 'bold');
    fprintf('  Visualization complete (showing %d/%d endpoints).\n', num_plot, num_total);
    catch ME
        warning('Could not create visualization plot: %s', ME.identifier, '%s', ME.message);
        fprintf('  Skipping visualization...\n');
    end
end

%% Write to CSV File
fprintf('Writing to CSV file: %s\n', outputFile);
fid = fopen(outputFile, 'w');
if fid == -1
    error('Could not open output file: %s', outputFile);
end

% Write header
fprintf(fid, 'coord_in_x,coord_in_y,coord_in_z,surf_in,coord_out_x,coord_out_y,coord_out_z,surf_out,parc_in,parc_out\n');

% Write data lines
num_streamlines = length(vtx_in);
for i = 1:num_streamlines
    % Format: coord_in(3), surf_in, coord_out(3), surf_out, parc_in, parc_out
    fprintf(fid, '%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%d,%d,%d\n', ...
        coords_in(i, 1), coords_in(i, 2), coords_in(i, 3), surf_in(i), ...
        coords_out(i, 1), coords_out(i, 2), coords_out(i, 3), surf_out(i), ...
        parc_in(i), parc_out(i));
    
    % Progress indicator for large files
    if mod(i, 10000) == 0
        fprintf('  Processed %d / %d streamlines...\n', i, num_streamlines);
    end
end

fclose(fid);
fprintf('Successfully wrote %d streamlines to %s\n', num_streamlines, outputFile);

% Print summary statistics
fprintf('\nSummary:\n');
fprintf('  Total streamlines: %d\n', num_streamlines);
fprintf('  Valid mappings (in): %d (%.1f%%)\n', sum(valid_in), 100*sum(valid_in)/num_streamlines);
fprintf('  Valid mappings (out): %d (%.1f%%)\n', sum(valid_out), 100*sum(valid_out)/num_streamlines);
fprintf('  Unique parcels (in): %d\n', length(unique(parc_in(valid_in))));
fprintf('  Unique parcels (out): %d\n', length(unique(parc_out(valid_out))));
fprintf('  Parcellation: %s\n', sbci_parc(parc_idx).atlas{1});
fprintf('  Grid file: %s\n', gridFile);

end

