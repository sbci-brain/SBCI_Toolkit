% MATLAB R2018a
%
% FUNCTION NAME:
%   plot_cortical
%
% DESCRIPTION:
%   Plot SBCI data on the surface given surface
%
% INPUT:
%   sbci_surf - (struct) A structure containing LH and RH meshes
%   sbci_map - (struct) A structure containing SBCI mapping information
%   data - (vector) A vector of data to plot according to the SBCI mapping
%   varargin - Optional arguments (legend=title to place over colorbar,
%      clim=limits for colorbar, cmap=colormap, bg=background color,
%      figid=id of the generated figure, logscale=display colorbar in log
%      scale, overlay_parc=parcellation structure for overlaying parcel
%      boundaries
%
% OUTPUT:
%   fig - (figure) handle to the resulting figure
%   ax - (axes) the 4 axes for the four viewpoints
%   cb - (colorbar) handle to the colorbar
%   Side effects: figure
%
% ASSUMPTIONS AND LIMITATIONS:
%   The surface, mapping, must come from the same run of the SBCI pipeline,
%   and the data must be a vector of the same length of the downsampled mesh 
%
function [fig, ax, cb] = plot_cortical(sbci_surf, sbci_map, data, varargin)

p = inputParser;
addParameter(p, 'legend', "", @ischar);
addParameter(p, 'clim', double([min(data) max(data)]), @isnumeric);
addParameter(p, 'cmap', 'jet', @ischar);
addParameter(p, 'bg', 'white', @ischar);
addParameter(p, 'figid', 1, @(n)validateattributes(n,{'numeric'},{'nonnegative'}));
addParameter(p, 'logscale', false, @islogical);
addParameter(p, 'overlay_parc',[], @isstruct)

% parse optional variables
parse(p, varargin{:});
params = p.Results;

% upsample data to the high resolution meshes
upsampled_data = zeros(sbci_map.shape(4), 1);

for i = 1:length(data)
    upsampled_data(sbci_map.map(1, sbci_map.map(2,:) == i)) = data(i);
end

% Upsample overlay labels if present 
if ~isempty(params.overlay_parc)
    upsampled_overlay = zeros(sbci_map.shape(4), 1); 
    for i = 1:length(params.overlay_parc.labels)
        upsampled_overlay(sbci_map.map(1, sbci_map.map(2,:) == i)) = params.overlay_parc.labels(i);
    end
end

% TODO: figure out nice values for these
h = 0.6;
w = 0.20;
m = 0.02;

fig = figure(params.figid);
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf = sbci_surf.lh_surf;
ax(1) = axes('position', [0.07 0.3 w h]);

%%%% plot LEFT hemisphere -90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

% Draw parcel boundaries
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay(1:sbci_map.shape(5)));
end

view(-90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading interp;

ax(2) = axes('position', [0.07+w+m 0.3 w h]);

%%%% plot LEFT hemisphere 90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

% Draw parcel boundaries
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay(1:sbci_map.shape(5)));
end
view(90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading interp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax(3) = axes('position', [0.07+2*(w+m) 0.3 w h]);
surf = sbci_surf.rh_surf;

%%%% plot RIGHT hemisphere -90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data((sbci_map.shape(5)+1):sbci_map.shape(4))), 'EdgeColor', 'none');

% Draw parcel boundaries
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay((sbci_map.shape(5)+1):sbci_map.shape(4)));
end
view(-90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading interp;

ax(4) = axes('position', [0.07+3*(w+m) 0.3 w h]);

%%%% plot RIGHT hemisphere 90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data((sbci_map.shape(5)+1):sbci_map.shape(4))), 'EdgeColor', 'none');

% Draw parcel boundaries
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay((sbci_map.shape(5)+1):sbci_map.shape(4)));
end
view(90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading interp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure all plots have the same (and valid) colour range
if params.clim(1) == params.clim(2)
    params.clim = params.clim(1) + [-1 0];
end

for i=1:length(ax) 
    set(ax(i), 'CLim', params.clim); 
    if params.logscale == true
        set(ax(i), 'ColorScale', 'log') 
    end
    colormap(ax(i), params.cmap); 
end 

% set a colour bar and place it at the bottom 
cb = colorbar('location', 'South');
set(cb, 'Position', [0.35 0.18 0.3 0.06]);
set(cb, 'XAxisLocation', 'bottom');

% set the title if given one
handle = get(cb, 'Title');
set(handle, 'String', params.legend);

% set the background colour and set the aspect ratio
set(gcf,'Color','white')
pos = get(gcf, 'Position');
set(gcf, 'Color', params.bg, 'InvertHardcopy', 'off', ...
    'Position', [pos(1), pos(2), pos(3), pos(3) / 2.5]);

end

% Function to draw parcel boundaries
function draw_boundary_lines(surf, labels_to_display)
    % Convert labels to per-triangle from per-vertex
    tri_labels = mode(labels_to_display(surf.tri), 2);
    
    % Determine edges and triangles
    edges = [surf.tri(:, [1, 2]); surf.tri(:, [2, 3]); surf.tri(:, [1, 3])];
    edge_triangles = [(1:size(surf.tri, 1))'; (1:size(surf.tri, 1))'; (1:size(surf.tri, 1))'];
    
    % Sort edges 
    [sorted_edges, sortIdx] = sortrows(sort(edges, 2));
    sorted_edge_triangles = edge_triangles(sortIdx);
    
    % Find pairs of triangles that share an edge
    diff_edges = diff(sorted_edges, 1, 1);
    edge_diffs = sum(diff_edges, 2) == 0;
    
    % Triangles pairs sharing the same edge are adjacent
    triangle_pairs = [sorted_edge_triangles([edge_diffs; false]), sorted_edge_triangles([false; edge_diffs])];
    
    % Filter pairs to those with different labels
    diff_labels = tri_labels(triangle_pairs(:,1)) ~= tri_labels(triangle_pairs(:,2));
    final_edges = sorted_edges([edge_diffs; false], :);
    final_edges = final_edges(diff_labels, :);
    
    % Draw lines for these boundary edges
    for i = 1:size(final_edges, 1)
        vtx_indices = final_edges(i, :);
        line(surf.vtx(1, vtx_indices), surf.vtx(2, vtx_indices), surf.vtx(3, vtx_indices), 'Color', 'k', 'LineWidth', 3);
    end
end