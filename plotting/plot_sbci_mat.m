% FUNCTION NAME:
%   plot_sbci_mat
%
% DESCRIPTION:
%   Plot a continuous SC or FC matrix according to a given parcellation
%
% INPUT:
%   data - (matrix) A PxP matrix of continuous connectivity data   
%   sbci_parc - (struct) A struct with parcellation output from SBCI
%   varargin - Optional arguments:
%       roi_mask - (vector) A vector of label IDs for ROIs to remove
%       triangular - (logical) If true, the FC and SC matrices are 
%           symmeterised before calculating SFC
%       roi_grid - (double) Alpha level for showing a grid lines for ROIs
%       clim - ([double]) Colorbar limits
%       legend - (string) Legend to put beside the colorbar
%
% OUTPUT:
%   fig - (figure) Handle to the generated figure
%   Side effects: figure
%
% ASSUMPTIONS AND LIMITATIONS:
%   Removes diagonals, assumes the SC and FC matrices are either
%   symmetric or triangular, and that the parcellation or SC, FC, 
%   matrices have not been rearranged in any way from SBCI output.
%
function [fig] = plot_sbci_mat(data, sbci_parc, varargin)

p = inputParser;
addParameter(p, 'triangular', false, @islogical);
addParameter(p, 'roi_mask', [], @isnumeric);
addParameter(p, 'roi_grid', 0, @isnumeric);
addParameter(p, 'clim', double([min(min(data)) max(max(data))]), @isnumeric);
addParameter(p, 'cmap', 'parula', @ischar);
addParameter(p, 'legend', "", @ischar);
addParameter(p, 'figid', 1, @(n)validateattributes(n,{'numeric'},{'nonnegative'}));

% parse optional variables
parse(p, varargin{:});
params = p.Results;
    
% symmeterise matrix
if (params.triangular == true)
    data = data + data';  
end
    
% remove diagonal elements
data = data - diag(diag(data));

lbl = sbci_parc.labels(sbci_parc.sorted_idx);
idx = sbci_parc.sorted_idx(~ismember(lbl, params.roi_mask));
data = data(idx, idx);

fig = figure(params.figid);
clf;

imagesc(data);

% display a grid of ROIs
[~,~,roi_labels] = unique(sbci_parc.labels(idx));
ticks = [0; find(diff(roi_labels)); size(data,1)];

% remove ticks
xticks([])
xticklabels([]);
yticks([])
yticklabels([]);

hold on

for i = 1:length(ticks)-1
    plot([ticks(i) ticks(i)],[ticks(1) ticks(end)], 'color', [0,0,0,params.roi_grid])
    plot([ticks(1) ticks(end)],[ticks(i) ticks(i)], 'color', [0,0,0,params.roi_grid])
end

hold off

axis square;
daspect([1,1,1]);

% ensure all plots have the same (and valid) colour range
if params.clim(1) == params.clim(2)
    params.clim = params.clim(1) + [-1 0];
end

set(gca, 'CLim', params.clim);
colormap(gca, params.cmap);

% set a colour bar
cb = colorbar();

% set the title if given one
handle = get(cb, 'Title');
set(handle, 'String', params.legend);

end