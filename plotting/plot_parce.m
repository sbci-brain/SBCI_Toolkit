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

p = inputParser;
addParameter(p, 'my_cols', []);
addParameter(p, 'overlay_parc',[], @isstruct)

% parse optional variables
parse(p, varargin{:});
params = p.Results;

% Get colors for each parcel
if isempty(params.my_cols)
    lab_color = distinguishable_colors(length(label));
else
    lab_color = params.my_cols; 
end
color_iter = 1; 

% upsample data to the high resolution meshes
upsampled_data = zeros(sbci_map.shape(4), 1);

for i = 1:length(label)
    upsampled_data(sbci_map.map(1, sbci_map.map(2,:) == i)) = label(i);
end

% Upsample overlay labels if present 
if ~isempty(params.overlay_parc)
    upsampled_overlay = zeros(sbci_map.shape(4), 1); 
    for i = 1:length(params.overlay_parc.labels)
        upsampled_overlay(sbci_map.map(1, sbci_map.map(2,:) == i)) = params.overlay_parc.labels(i);
    end
end

%% %%%%%%%%%%left surface plot 1 %%%%%%%%%%%% %%
surf = sbci_surf.lh_surf;
subplot(2,2,1)
ax1 = gca; 
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),'EdgeColor',...
    'none','Facecolor', 'k');
title('Left Hemisphere - Front')
daspect([1 1 1]); 
axis tight; 
view(90,23)

labels_to_display = upsampled_data(1:sbci_map.shape(5));
unique_labs = unique(labels_to_display);

for i = 1:length(unique_labs)
    hold on

    %find tri corresponding to unique_labs(i)
    vertex_lab = find(labels_to_display == unique_labs(i));  
    for j=1:3
        is_member(j,:) = ismember(surf.tri(:,j),vertex_lab);
    end
    surf_idx = sum(is_member)>1; 

    trimesh(surf.tri(surf_idx,:), surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
       'FaceColor', lab_color(color_iter,:), 'EdgeColor',lab_color(color_iter,:),'FaceAlpha',0.8);
    color_iter = color_iter + 1; 
end
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay(1:sbci_map.shape(5)));
end
axis off;

%% %%%%%%%%%%left surface plot 2 %%%%%%%%%%%% %%
subplot(2,2,3)
ax2 = copyobj(ax1,gcf); 
set(ax2, 'Position', get(subplot(2,2,3), 'Position'));
title(ax2, 'Left Hemisphere - Back')
view(ax2, 270, 0)
axis off;

%% %%%%%%%%%%right surface plot 1 %%%%%%%%%%%% %% 

surf = sbci_surf.rh_surf;

subplot(2,2,2)
ax3 = gca; 
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),'EdgeColor',...
    'none','Facecolor', 'k');
daspect([1 1 1]); 
axis tight; 
view(-90,23)
title('Right Hemisphere - Front')

labels_to_display = upsampled_data((sbci_map.shape(5)+1):sbci_map.shape(4));
unique_labs = unique(labels_to_display);

for i = 1:length(unique_labs)
    hold on
    
    %find tri corresponding to unique_labs(i)
    vertex_lab = find(labels_to_display == unique_labs(i));
    
    for j=1:3
        is_member(j,:) = ismember(surf.tri(:,j),vertex_lab);
    end
    
    surf_idx = sum(is_member)>1;   
    trimesh(surf.tri(surf_idx,:), surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
       'FaceColor', lab_color(color_iter,:), 'EdgeColor',lab_color(color_iter,:),'FaceAlpha',0.8);
    color_iter = color_iter + 1; 
    
end
if ~isempty(params.overlay_parc)
    draw_boundary_lines(surf, upsampled_overlay((sbci_map.shape(5)+1):sbci_map.shape(4)));
end
axis off;

%% %%%%%%%%%%right surface plot 2 %%%%%%%%%%%% %% 

subplot(2,2,4)
ax4 = copyobj(ax3,gcf); 
set(ax4, 'Position', get(subplot(2,2,4), 'Position'));
title(ax4, 'Right Hemisphere - Back')
view(ax4, 90, 0)
axis off;

set(gcf, 'Color', 'w');
result = 1;
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
