% MATLAB R2018a
%
% FUNCTION NAME:
%   load_sbci_data
%
% DESCRIPTION:
%   Loads required data for analysis of SBCI data
%
% INPUT:
%   sbci_location - (string) location of the SBCI_AVG output folder
%   resolution - (float) value between 0-1 corresponding the the resolution
%       of SBCI mappings required.
%
% OUTPUT:
%   sbci_parc - (struct array) structures containing all the processed
%       parcellation informatiob at the given resolution
%   sbci_mapping - (struct) structure containing all the SBCI mapping
%       information from high resolution to the resolution given
%   adjacency - (mat) and adjacency matrix for the given resolution
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   Assumes the SBCI pipeline ran successfully 
%
function [sbci_parc, sbci_mapping, adjacency] = load_sbci_data(sbci_location, resolution)

% load the mapping and adjacency matrix at the given resolution
sbci_mapping = load(sprintf('%s/mapping_avg_%.2f.mat', sbci_location, resolution));
adjacency = load(sprintf('%s/adjacency_%.2f.mat', sbci_location, resolution));

% find all the parcellation files
files = dir(sprintf('%s/*_avg_roi_%.2f.mat', sbci_location, resolution));

for i=1:length(files)
    parc = load(files(i).name);
    parc_name = strsplit(files(i).name, '_avg_roi');
    
    % put parcellations in a vector of structs
    sbci_parc(i).atlas = parc_name(1);
    sbci_parc(i).sorted_idx = parc.sorted_idx;
    sbci_parc(i).labels = parc.labels;
    sbci_parc(i).names = parc.names;
end

end
