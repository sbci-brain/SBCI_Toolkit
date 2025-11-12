% MATLAB R2018a
%
% FUNCTION NAME:
%   plot_parce_with_endpoints
%   
% DESCRIPTION:
%   Visualize parcellation with transparency and overlay streamline endpoints
%
% INPUT:
%   sbci_surf - (struct) A structure containing LH and RH meshes
%   sbci_map - (struct) A structure containing SBCI mapping information
%   labels - (vector) A vector with parcellation labels according to SBCI
%       format, corresponding to the connectome cols or rows
%   coords_in - (matrix) [Nx3] Cartesian coordinates of streamline start points
%   coords_out - (matrix) [Nx3] Cartesian coordinates of streamline end points
%   surf_in - (vector) [Nx1] Hemisphere indicator for start points (0=left, 1=right)
%   surf_out - (vector) [Nx1] Hemisphere indicator for end points (0=left, 1=right)
%   varargin - Optional arguments:
%       my_cols - (matrix) A Kx3 vector of colors for each parcel
%       parcel_alpha - (scalar) Transparency level for parcels (0-1, default 0.5)
%       endpoint_size - (scalar) Size of endpoint markers (default 20)
%       show_lines - (logical) Whether to draw lines connecting endpoint pairs (default false)
%
% SIDE-EFFECT: plot of the parcellation with endpoints
% OUTPUT:
%   1
%
% ASSUMPTIONS AND LIMITATIONS:
%   The data and mapping must come from the same run of the SBCI pipeline.

function result = plot_parce_with_endpoints(sbci_surf, sbci_map, label, ...
    coords_in, coords_out, surf_in, surf_out, varargin)

% Parse input 
p = inputParser;
addParameter(p, 'my_cols', []);
addParameter(p, 'parcel_alpha', 0.5, @isnumeric);
addParameter(p, 'endpoint_size', 20, @isnumeric);
addParameter(p, 'show_lines', false, @islogical);
addParameter(p, 'selected_parcels', [], @isnumeric);
parse(p, varargin{:});
params = p.Results;

% If user didn't provide a colormap, create a default
if isempty(params.my_cols)
    % Get all unique labels
    unique_labels = unique(label);
    unique_labels = unique_labels(unique_labels ~= 0);  % Remove 0 if present
    if isempty(unique_labels)
        unique_labels = [1];  % Fallback
    end
    
    % Create colormap for all parcels
    lab_color = distinguishable_colors(numel(unique_labels));
    
    % If selected_parcels is provided, we'll handle gray coloring in the mapping
    % Otherwise, use normal colors for all
else
    lab_color = params.my_cols;
end

% Upsample labels to high-res vertices 
total_verts = sbci_map.shape(4);
upsampled_data = zeros(total_verts, 1);

vertex_indices = sbci_map.map(1, :);
label_indices  = sbci_map.map(2, :);

upsampled_data(vertex_indices) = label(label_indices);

% Split LH / RH data
n_lh = sbci_map.shape(5);          
lh_labels = upsampled_data(1:n_lh);
rh_labels = upsampled_data(n_lh+1:end);

% Set figure background to white
figure('Color', 'w');


%% ========== LEFT HEMISPHERE ==========
% LH surface 
surf_lh = sbci_surf.lh_surf;  % .vtx is 3 x N, .tri is M x 3

% Labels for each triangle
tri_labels_lh = mode(lh_labels(surf_lh.tri), 2);
% Map labels to sequential indices for colormap
% If selected_parcels is provided, unselected parcels will be gray
if ~isempty(params.selected_parcels)
    % Separate selected and unselected parcels
    selected_mask_lh = ismember(tri_labels_lh, params.selected_parcels);
    unselected_mask_lh = ~selected_mask_lh & (tri_labels_lh ~= 0);
    
    % Get unique selected labels
    unique_selected = unique(tri_labels_lh(selected_mask_lh));
    if isempty(unique_selected)
        unique_selected = [1];  % Fallback
    end
    
    % Create mapping: unselected -> 1 (gray), selected -> 2, 3, 4, ...
    [~, label_idx] = ismember(tri_labels_lh, unique_selected);
    tri_labels_lh(selected_mask_lh) = label_idx(selected_mask_lh) + 1;  % Selected parcels start at index 2
    tri_labels_lh(unselected_mask_lh) = 1;  % Unselected parcels -> gray (index 1)
    tri_labels_lh(tri_labels_lh == 0) = 1;  % Zero labels -> gray
    
    % Update colormap: gray at index 1, then colors for selected parcels
    lab_color = [0.7, 0.7, 0.7; distinguishable_colors(numel(unique_selected))];
else
    % No selection - use normal mapping
    unique_labels_nonzero = unique(tri_labels_lh(tri_labels_lh ~= 0));
    if isempty(unique_labels_nonzero)
        unique_labels_nonzero = [1];
    end
    [~, label_idx] = ismember(tri_labels_lh, unique_labels_nonzero);
    tri_labels_lh(tri_labels_lh ~= 0) = label_idx(tri_labels_lh ~= 0) + 1;
    tri_labels_lh(tri_labels_lh == 0) = 1;
    lab_color = [0.7, 0.7, 0.7; lab_color];  % Gray at index 1
end

% Plot LH in subplot(2,2,1)
subplot(2,2,1);
hold on;
title('Left Hemisphere - Front');

% Create patch with transparency
% All parcels use the same alpha (no transparency for unselected)
face_alpha_lh = params.parcel_alpha * ones(size(tri_labels_lh));

h_patch_lh1 = patch( ...
    'Vertices', surf_lh.vtx', ...
    'Faces',    surf_lh.tri, ...
    'FaceVertexCData', tri_labels_lh, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'CDataMapping', 'direct', ...
    'FaceVertexAlphaData', face_alpha_lh, ...
    'FaceAlpha', 'flat');

daspect([1 1 1]);
axis off;
axis tight;
view(90,23);  % lateral view
colormap(lab_color);  % use your desired colormap

% Overlay endpoints on left hemisphere
left_in_idx = find(surf_in == 0);
left_out_idx = find(surf_out == 0);

% Plot start points on left hemisphere
if ~isempty(left_in_idx) && ~isempty(coords_in)
    scatter3(coords_in(left_in_idx,1), coords_in(left_in_idx,2), coords_in(left_in_idx,3), ...
        params.endpoint_size, [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% Plot end points on left hemisphere
if ~isempty(left_out_idx) && ~isempty(coords_out)
    scatter3(coords_out(left_out_idx,1), coords_out(left_out_idx,2), coords_out(left_out_idx,3), ...
        params.endpoint_size, [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% Draw lines connecting endpoint pairs if requested
if params.show_lines
    % Find pairs where both endpoints are in left hemisphere
    left_pairs = (surf_in == 0) & (surf_out == 0);
    if any(left_pairs)
        for i = find(left_pairs)'
            line([coords_in(i,1), coords_out(i,1)], ...
                 [coords_in(i,2), coords_out(i,2)], ...
                 [coords_in(i,3), coords_out(i,3)], ...
                 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        end
    end
end

%%% Left Hemisphere Plot 2 (Back View) %%%
ax1 = gca;  
subplot(2,2,3);
ax2 = copyobj(ax1, gcf); 
set(ax2, 'Position', get(subplot(2,2,3), 'Position'));
title(ax2, 'Left Hemisphere - Back');
view(ax2, [270, 0]);   % back view
axis off;

%% ========== RIGHT HEMISPHERE ==========
surf_rh = sbci_surf.rh_surf;  % .vtx is 3 x N, .tri is M x 3

tri_labels_rh = mode(rh_labels(surf_rh.tri), 2);
% Map labels to sequential indices for colormap (same logic as LH)
if ~isempty(params.selected_parcels)
    % Separate selected and unselected parcels
    selected_mask_rh = ismember(tri_labels_rh, params.selected_parcels);
    unselected_mask_rh = ~selected_mask_rh & (tri_labels_rh ~= 0);
    
    % Get unique selected labels
    unique_selected = unique(tri_labels_rh(selected_mask_rh));
    if isempty(unique_selected)
        unique_selected = [1];  % Fallback
    end
    
    % Create mapping: unselected -> 1 (gray), selected -> 2, 3, 4, ...
    [~, label_idx] = ismember(tri_labels_rh, unique_selected);
    tri_labels_rh(selected_mask_rh) = label_idx(selected_mask_rh) + 1;  % Selected parcels start at index 2
    tri_labels_rh(unselected_mask_rh) = 1;  % Unselected parcels -> gray (index 1)
    tri_labels_rh(tri_labels_rh == 0) = 1;  % Zero labels -> gray
else
    % No selection - use normal mapping
    unique_labels_nonzero = unique(tri_labels_rh(tri_labels_rh ~= 0));
    if isempty(unique_labels_nonzero)
        unique_labels_nonzero = [1];
    end
    [~, label_idx] = ismember(tri_labels_rh, unique_labels_nonzero);
    tri_labels_rh(tri_labels_rh ~= 0) = label_idx(tri_labels_rh ~= 0) + 1;
    tri_labels_rh(tri_labels_rh == 0) = 1;
end

subplot(2,2,2);
hold on;
title('Right Hemisphere - Front');

% Create patch with transparency
% All parcels use the same alpha (no transparency for unselected)
face_alpha_rh = params.parcel_alpha * ones(size(tri_labels_rh));

h_patch_rh1 = patch( ...
    'Vertices', surf_rh.vtx', ...
    'Faces',    surf_rh.tri, ...
    'FaceVertexCData', tri_labels_rh, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'CDataMapping', 'direct', ...
    'FaceVertexAlphaData', face_alpha_rh, ...
    'FaceAlpha', 'flat');

daspect([1 1 1]);
axis off;
axis tight;
view(-90,23);  % lateral view
colormap(lab_color);

% Overlay endpoints on right hemisphere
right_in_idx = find(surf_in == 1);
right_out_idx = find(surf_out == 1);

% Plot start points on right hemisphere
if ~isempty(right_in_idx) && ~isempty(coords_in)
    scatter3(coords_in(right_in_idx,1), coords_in(right_in_idx,2), coords_in(right_in_idx,3), ...
        params.endpoint_size, [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% Plot end points on right hemisphere
if ~isempty(right_out_idx) && ~isempty(coords_out)
    scatter3(coords_out(right_out_idx,1), coords_out(right_out_idx,2), coords_out(right_out_idx,3), ...
        params.endpoint_size, [0 0 0], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% Draw lines connecting endpoint pairs if requested
if params.show_lines
    % Find pairs where both endpoints are in right hemisphere
    right_pairs = (surf_in == 1) & (surf_out == 1);
    if any(right_pairs)
        for i = find(right_pairs)'
            line([coords_in(i,1), coords_out(i,1)], ...
                 [coords_in(i,2), coords_out(i,2)], ...
                 [coords_in(i,3), coords_out(i,3)], ...
                 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        end
    end
end

%%% Right Hemisphere Plot 2 (Back View) %%%
ax3 = gca;
subplot(2,2,4);
ax4 = copyobj(ax3, gcf);
set(ax4, 'Position', get(subplot(2,2,4), 'Position'));
title(ax4, 'Right Hemisphere - Back');
view(ax4, [90, 0]);
axis off;

% No legend needed - endpoints are all the same color

result = 1;
end

