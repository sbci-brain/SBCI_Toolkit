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
%      clim=limits for colorbar, cmap=colormap, bg=background color
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

% parse optional variables
parse(p, varargin{:});
params = p.Results;

% upsample data to the high resolution meshes
upsampled_data = zeros(sbci_map.shape(4), 1);

for i = 1:length(data)
    upsampled_data(sbci_map.map(1, sbci_map.map(2,:) == i)) = data(i);
end

% TODO: figure out nice values for these
h = 0.6;
w = 0.20;
m = 0.02;

fig = figure();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf = sbci_surf.lh_surf;
ax(1) = axes('position', [0.07 0.3 w h]);

%%%% plot LEFT hemisphere -90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

view(-90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading flat;

ax(2) = axes('position', [0.07+w+m 0.3 w h]);

%%%% plot LEFT hemisphere 90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

view(90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading flat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax(3) = axes('position', [0.07+2*(w+m) 0.3 w h]);
surf = sbci_surf.rh_surf;

%%%% plot RIGHT hemisphere -90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

view(-90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading flat;

ax(4) = axes('position', [0.07+3*(w+m) 0.3 w h]);

%%%% plot RIGHT hemisphere 90 %%%%
trisurf(surf.tri, surf.vtx(1,:), surf.vtx(2,:), surf.vtx(3,:),...
   double(upsampled_data(1:sbci_map.shape(5))), 'EdgeColor', 'none');

view(90,0)

% options to make the plot look pretty
daspect([1 1 1]); axis tight; camlight; axis vis3d off;
lighting gouraud; material dull; shading flat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ensure all plots have the same (and valid) colour range
if params.clim(1) == params.clim(2)
    params.clim = params.clim(1) + [-1 0];
end

for i=1:length(ax)
    set(ax(i), 'CLim', params.clim);
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
whitebg(gcf, params.bg);
pos = get(gcf, 'Position');
set(gcf, 'Color', params.bg, 'InvertHardcopy', 'off', ...
    'Position', [pos(1), pos(2), pos(3), pos(3) / 2.5]);

end