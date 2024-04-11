% FUNCTION NAME:
%   plot_sfc
%
% DESCRIPTION:
%   Plots continuous FC and SC matrices as per the specified parcellation.
%
% INPUT:
%   fc - (matrix) A PxP matrix of continuous functional connectivity data.
%   sc - (matrix) A PxP matrix of continuous structural connectivity data.
%   sbci_parc - (struct) A struct with parcellation output from SBCI.
%   atlas_num - The number of the atlas to use for parcellation, which indexes into sbci_parc.
%   roi_mask_num - Regions that not interest to and want to removed from
%   the visualization.
% OUTPUT:
%   Generates two figures:
%   The first figure represents the continuous FC matrix by the given resolution.
%   The second figure represents the continuous SC matrix by the given resolution.
%   Side effects: Two figures (plot windows) are opened, displaying the respective connectivity matrices.
% ASSUMPTIONS AND LIMITATIONS:
%   None


function plot_connectivity(fc, sc, sbci_parc, atlas_index, roi_exclusion_index)
    % convert into full matrices (for plotting)
    fc = fc + fc' - 2*diag(diag(fc));
    sc = sc + sc' - 2*diag(diag(sc));
    sc = sc/sum(sum(sc));

    plot_sbci_mat(fc, sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 1, 'clim', [-0.1, 0.35]);
    title(['Continuous FC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');

    plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 2, 'clim', [0, 3.5]);
    title(['Continuous SC (' sbci_parc(atlas_index).atlas{1} ')'], 'Interpreter', 'none');
end