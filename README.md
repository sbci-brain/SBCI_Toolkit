# SBCI_Toolkit
MATLAB scripts for visualisation and analysis of SBCI data

## Basic usage

```
% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('`LOCATION_OF_SBCI_AVG_DATA`', `RESOLUTION`);
 
% example to plot data on the inflated cortical surface
sbci_surf = load_sbci_surfaces('`LOCATION_OF_SBCI_AVG_DATA`');
plot_cortical(sbci_surf.inflated, sbci_mapping, `DATA`);

% example to plot SC or FC matrix with a given parcellation
plot_sbci_mat(`DATA`, sbci_parc(0), 'roi_mask', [1,35]);
```

For the processed HCP data, `RESLOLUTION`=0.94

`sbci_parc` is an array of structures, each one corresponding to a different parcellation. Here, I'm assuming that sbci_parc(0) corresponds to the Desikan atlas, but 
you can check the names by looking at `sbci_parc(0).atlas`, for example. As for the `'roi_mask'`, this will remove the selected regions from the matrix plot, and label
names can also be checked by looking at `sbci_parc(0).names`.