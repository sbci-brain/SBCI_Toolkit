% FUNCTION NAME:
%   plot_value_cortically
%
% DESCRIPTION:
%   This function loads numerical values from a specified text file and visualizes 
%   these values on an inflated cortical surface model provided by the user.
%
% INPUT:
%   sbci_surf - (struct) A struct containing cortical surface information,
%   including an inflated cortical surface model.
%   sbci_mapping - (mapping) Mapping information used for the parcellation of the FC and SC matrices.
%   txt_file - (string) The path to the text file containing the numerical 
%   values to be visualized. The file should contain a list or matrix of values 
%   that correspond to specific regions or vertices on the cortical surface.
% OUTPUT:
%   Generates a figure displaying the loaded values mapped onto the cortical surface. 
% ASSUMPTIONS AND LIMITATIONS:
%   None

function plot_value_cortically(sbci_surf, sbci_mapping, txt_file)
   txt = load(txt_file);
   plot_cortical(sbci_surf.inflated, sbci_mapping, txt);
end