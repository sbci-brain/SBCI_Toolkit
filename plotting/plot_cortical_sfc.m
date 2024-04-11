% FUNCTION NAME:
%   plot_cortical_sfc
%
% DESCRIPTION:
%   The function calculates global SFC by considering the entire brain connectivity and local SFC based on specific parcellation regions, then visualizes these metrics on cortical surfaces.
%
% INPUT:
%   fc - (matrix) A PxP matrix of continuous functional connectivity data.
%   sc - (matrix) A PxP matrix of continuous structural connectivity data.
%   sbci_surf - (struct) A struct containing cortical surface information, including an inflated cortical surface model.
%   sbci_mapping - (mapping) Mapping information used for the parcellation of the FC and SC matrices.
%   sbci_parc - (struct) A struct with parcellation output from SBCI.
%   atlas_num - The number of the atlas to use for parcellation, which indexes into sbci_parc.
% OUTPUT:
%   Generates two figures:
%   The first figure displays the global SFC mapped onto the cortical surface.
%   The second figure displays the local SFC (based on specified parcellation) on the cortical surface.
%   Side effects: Two figures (plot windows) are opened, displaying the cortical surface.
% ASSUMPTIONS AND LIMITATIONS:
%   None

function plot_cortical_sfc(sfc_gbl, sfc_loc, sbci_surf, sbci_mapping, sbci_parc, atlas_num)

    sfc_gbl(isnan(sfc_gbl)) = 0;
    sfc_loc(isnan(sfc_loc)) = 0;

    plot_cortical(sbci_surf.inflated, sbci_mapping, sfc_gbl, 'figid', 1, ...
      'legend', 'SFC global', 'clim', [0,1])

    plot_cortical(sbci_surf.inflated, sbci_mapping, sfc_loc, 'figid', 2, ...
      'legend', ['SFC local (' sbci_parc(atlas_num).atlas{1} ')'], 'clim', [0,1])
end