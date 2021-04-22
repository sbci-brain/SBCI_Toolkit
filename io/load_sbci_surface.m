% FUNCTION NAME:
%   load_sbci_surfaces
%
% DESCRIPTION:
%   Loads required SBCI surfaces for analysis (mappings, surfaces, etc.)
%
% INPUT:
%   sbci_location - (string) location of the SBCI_AVG output folder
%
% OUTPUT:
%   sbci_surf - (struct) a structure with all the surfaces (inflated, etc.)
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   Assumes the SBCI pipeline ran successfully 
%
function [sbci_surf] = load_sbci_surface(sbci_location)

sbci_surf.inflated.lh_surf = read_vtk(sprintf('%s/lh_inflated_avg_lps.vtk', sbci_location));
sbci_surf.inflated.rh_surf = read_vtk(sprintf('%s/rh_inflated_avg_lps.vtk', sbci_location));

sbci_surf.white.lh_surf = read_vtk(sprintf('%s/lh_white_avg_lps.vtk', sbci_location));
sbci_surf.white.rh_surf = read_vtk(sprintf('%s/rh_white_avg_lps.vtk', sbci_location));

sbci_surf.sphere.lh_surf = read_vtk(sprintf('%s/lh_sphere_avg_lps.vtk', sbci_location));
sbci_surf.sphere.rh_surf = read_vtk(sprintf('%s/rh_sphere_avg_lps.vtk', sbci_location));

end
