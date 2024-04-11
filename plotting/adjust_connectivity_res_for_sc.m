% FUNCTION NAME:
%   adjust_sfc_res
%
% DESCRIPTION:
%   visualizes the adjusted SC matrices using parcellation information
%
% INPUT:
%   sc - (matrix) A PxP matrix of continuous structural connectivity data.
%   sbci_parc - (struct) A struct with parcellation output from SBCI
%   atlas_num - The number of the atlas to use for parcellation, which indexes into sbci_parc.
%   sbci_mapping - (struct) A structure containing SBCI mapping information
%   roi_mask_num - (vector) A vector of label IDs for ROIs to remove
% OUTPUT:
%   A figure displays the discrete FC matrix with given visualization parameters.
% ASSUMPTIONS AND LIMITATIONS:
%   None

function adjust_connectivity_res_for_sc(sc, sbci_parc, atlas_num, sbci_mapping, roi_mask_num)

    % Parcellate the FC/SC matrix based on given atlas
    dct_sc = parcellate_sc(sc, sbci_parc(atlas_num), sbci_mapping, 'roi_mask', roi_mask_num);
    
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