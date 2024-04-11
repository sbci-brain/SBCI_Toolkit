% FUNCTION NAME:
%   adjust_connectivity_res_for_fc
%
% DESCRIPTION:
%   visualizes the adjusted FC and SC matrices using parcellation information
%
% INPUT:
%   fc - (matrix) A PxP matrix of continuous functional connectivity data.
%   sbci_parc - (struct) A struct with parcellation output from SBCI
%   atlas_num - The number of the atlas to use for parcellation, which indexes into sbci_parc.
%   sbci_mapping - (struct) A structure containing SBCI mapping information
%   roi_mask_num - (vector) A vector of label IDs for ROIs to remove
% OUTPUT:
%   A figure displays the discrete FC matrix with given visualization parameters.
% ASSUMPTIONS AND LIMITATIONS:
%   None
function adjust_connectivity_res_for_fc(fc, sbci_parc, atlas_num, sbci_mapping, roi_mask_num)
    % Parcellate the FC/SC matrix based on given atlas
    dct_fc = parcellate_fc(fc, sbci_parc(atlas_num), sbci_mapping, 'roi_mask', roi_mask_num);
    
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
end