% FUNCTION NAME:
%   adjust_sfc_res
%
% DESCRIPTION:
%   visualizes the adjusted FC and SC matrices using parcellation information
%
% INPUT:
%   fc - (matrix) A PxP matrix of continuous functional connectivity data.
%   sc - (matrix) A PxP matrix of continuous structural connectivity data.
%   sbci_parc - (struct) A struct with parcellation output from SBCI.
%   atlas_num - The number of the atlas to use for parcellation, which indexes into sbci_parc.
%   sbci_mapping - (mapping) Mapping information used for the parcellation of the FC and SC matrices.
%   roi_mask_num - Regions that not interest to and want to removed from
%   the visualization.
% OUTPUT:
%   Generates two figures:
%   The first figure displays the discrete FC matrix with given visualization parameters.
%   The second figure displays the discrete SC matrix with given visualization parameters.
%   Side effects: Two figures (plot windows) are opened, displaying the respective connectivity matrices.
% ASSUMPTIONS AND LIMITATIONS:
%   None

function adjust_connectivity_res(fc, sc, sbci_parc, atlas_num, sbci_mapping, roi_mask_num)
    % convert into full matrices (for plotting)
    fc = fc + fc' - 2*diag(diag(fc));
    sc = sc + sc' - 2*diag(diag(sc));
    sc = sc/sum(sum(sc));

    % Parcellate the FC/SC matrix based on given atlas
    dct_fc = parcellate_fc(fc, sbci_parc(atlas_num), sbci_mapping, 'roi_mask', roi_mask_num);
    dct_sc = parcellate_sc(sc, sbci_parc(atlas_num), sbci_mapping, 'roi_mask', roi_mask_num);
    
    % Plot Discrete FC
    figure(3);
    imagesc(dct_fc);
    
    % remove ticks
    xticks([]); xticklabels([]);
    yticks([]); yticklabels([]);
    
    % colours and aspect ratio
    axis square;
    daspect([1,1,1]);
    set(gca, 'CLim', [-0.1, 0.35]);
    colorbar();
    
    title(['Discrete FC (' sbci_parc(atlas_num).atlas{1} ')'], 'Interpreter', 'none')
    
    % Plot Discrete SC
    figure(4);
    imagesc(log((10^7*dct_sc) + 1));

    % remove ticks
    xticks([]); xticklabels([]);
    yticks([]); yticklabels([]);
    
    % colours and aspect ratio
    axis square;
    daspect([1,1,1]);
    set(gca, 'CLim', [0, 3.5]);
    colorbar();
    
    title(['Discrete SC (' sbci_parc(atlas_num).atlas{1} ')'], 'Interpreter', 'none')
end