% MATLAB R2018a
%
% FUNCTION NAME:
%   plot_parce
%   
% DESCRIPTION:
%   visualize parcellation
%
% INPUT:
%   sbci_surf - sbci_surf - (struct) A structure containing LH and RH meshes
%   sbci_map - (struct) A structure containing SBCI mapping information
%   labels - (vector) A vector with parcellation labels according to SBCI
%   format, corresponding to the connectome cols or rows
%   varargin - Optional arguments:
%       my_cols - (vector) A Kx3 vector of colors for each parcel
%       overlay_parc - (struct) A parcellation structure for overlaying
%       parcel boundaries 
%
% SIDE-EFFECT: plot of the parcellation
% OUTPUT:
%   1
%
% ASSUMPTIONS AND LIMITATIONS:
%   The data and mapping must come from the same run of the SBCI pipeline.
%

function result = plot_parce(sbci_surf, sbci_map, label,varargin)

% Parse input 
p = inputParser;
addParameter(p, 'my_cols', []);
addParameter(p, 'overlay_parc',[], @isstruct)
parse(p, varargin{:});
params = p.Results;

% If user didn't provide a colormap, create a default
if isempty(params.my_cols)
    unique_labels = unique(label);
    lab_color = distinguishable_colors(numel(unique_labels));
else
    lab_color = params.my_cols;
end

% Upsample labels to high-res vertices 
total_verts = sbci_map.shape(4);
upsampled_data = zeros(total_verts, 1);

vertex_indices = sbci_map.map(1, :);
label_indices  = sbci_map.map(2, :);

upsampled_data(vertex_indices) = label(label_indices);

% If there's an overlay parcellation
if ~isempty(params.overlay_parc)
    upsampled_overlay = zeros(total_verts, 1);
    upsampled_overlay(vertex_indices) = params.overlay_parc.labels(label_indices);
end

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
[~, ~, ic] = unique(tri_labels_lh, 'sorted');  
tri_labels_lh = ic;

% Plot LH in subplot(2,2,1)
subplot(2,2,1);
hold on;
title('Left Hemisphere - Front');

% single patch call
patch( ...
    'Vertices', surf_lh.vtx', ...
    'Faces',    surf_lh.tri, ...
    'FaceVertexCData', tri_labels_lh, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'CDataMapping', 'direct');

daspect([1 1 1]);
axis off;
axis tight;
view(90,23);  % lateral view
colormap(lab_color);  % use your desired colormap

% Overlay boundary lines if requested
if ~isempty(params.overlay_parc)
    lh_overlay = upsampled_overlay(1:n_lh);
    draw_boundary_lines(surf_lh, lh_overlay);
end

%%% Left Hemisphere Plot 2 %%%
ax1 = gca;  
subplot(2,2,3);
ax2 = copyobj(ax1, gcf); 
set(ax2, 'Position', get(subplot(2,2,3), 'Position'));
title(ax2, 'Left Hemisphere - Back');
view(ax2, [270, 0]);   % e.g., 270 azimuth, 0 elevation
axis off;

%% ========== RIGHT HEMISPHERE ==========
surf_rh = sbci_surf.rh_surf;  % .vtx is 3 x N, .tri is M x 3

tri_labels_rh = mode(rh_labels(surf_rh.tri), 2);
[~, ~, ic] = unique(tri_labels_rh, 'sorted');  
tri_labels_rh = ic;

subplot(2,2,2);
hold on;
title('Right Hemisphere - Front');

patch( ...
    'Vertices', surf_rh.vtx', ...
    'Faces',    surf_rh.tri, ...
    'FaceVertexCData', tri_labels_rh, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'CDataMapping', 'direct' ...
    );

daspect([1 1 1]);
axis off;
axis tight;
view(-90,23);  % lateral view
colormap(lab_color);

if ~isempty(params.overlay_parc)
    rh_overlay = upsampled_overlay(n_lh+1:end);
    draw_boundary_lines(surf_rh, rh_overlay);
end

result = 1;


%%% Right Hemisphere Plot 2 %%%
ax3 = gca;
subplot(2,2,4);
ax4 = copyobj(ax3, gcf);
set(ax4, 'Position', get(subplot(2,2,4), 'Position'));
title(ax4, 'Right Hemisphere - Back');
view(ax4, [90, 0]);
axis off
result = 1;
end

%%---------------------------------------------------------------------%%
% Function to draw boundaries between adjacent triangles with different

function draw_boundary_lines(surf, labels_of_vertices)

% Ensure we do not overwrite the patch
hold on;

% 1) Per-triangle label
tri_labels = mode(labels_of_vertices(surf.tri), 2);

% 2) Gather edges for each triangle
edges = [ ...
    surf.tri(:, [1,2]);
    surf.tri(:, [2,3]);
    surf.tri(:, [1,3]);
    ];
edge_triangles = [ ...
    (1:size(surf.tri, 1))';
    (1:size(surf.tri, 1))';
    (1:size(surf.tri, 1))';
    ];

% 3) Sort edges so duplicates line up
[sorted_edges, sortIdx] = sortrows(sort(edges, 2));
sorted_edge_triangles   = edge_triangles(sortIdx);

% 4) Identify shared edges (duplicates)
duplicate_flags = all(diff(sorted_edges) == 0, 2);

% Triangles that share the same edge
tri_pair_1 = sorted_edge_triangles([duplicate_flags; false]);
tri_pair_2 = sorted_edge_triangles([false; duplicate_flags]);
shared_edges = sorted_edges([duplicate_flags; false], :);

% 5) Keep only edges where the two adjacent triangles have different labels
is_diff_label = tri_labels(tri_pair_1) ~= tri_labels(tri_pair_2);
final_edges   = shared_edges(is_diff_label, :);

% 6) Build a big list of line endpoints
X = [surf.vtx(1, final_edges(:,1)); surf.vtx(1, final_edges(:,2))];
Y = [surf.vtx(2, final_edges(:,1)); surf.vtx(2, final_edges(:,2))];
Z = [surf.vtx(3, final_edges(:,1)); surf.vtx(3, final_edges(:,2))];

X = reshape(X, 2, []);
Y = reshape(Y, 2, []);
Z = reshape(Z, 2, []);

% Plot them all in one shot
line(X, Y, Z, 'Color', 'k', 'LineWidth', 2);
end
