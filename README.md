# SBCI_Toolkit
MATLAB scripts for visualisation and analysis of SBCI data

## Basic usage

```
% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('LOCATION_OF_SBCI_AVG_DATA', RESOLUTION);
 
% example to plot data on the inflated cortical surface
sbci_surf = load_sbci_surfaces('LOCATION_OF_SBCI_AVG_DATA');
plot_cortical(sbci_surf.inflated, sbci_mapping, DATA);
```

![SFC_example](https://user-images.githubusercontent.com/82663099/115779783-63c91600-a386-11eb-8b2d-e06d4384f562.png)

```
% example to plot SC or FC matrix with a given parcellation
plot_sbci_mat(DATA, sbci_parc(3), 'roi_mask', [1,36]);
```

![SC_FC_example](https://user-images.githubusercontent.com/82663099/115779808-6af02400-a386-11eb-8c0a-cf1ce5eff6e6.png)

For the processed HCP data, `RESLOLUTION`=0.94

`sbci_parc` is an array of structures, each one corresponding to a different parcellation. Here, I'm assuming that sbci_parc(3) corresponds to the Desikan atlas, but 
you can check the names by looking at `sbci_parc(3).atlas`, for example. As for the `'roi_mask'`, this will remove the selected regions from the matrix plot, and label
names can also be checked by looking at `sbci_parc(3).names`.

## Changing Resolution

SBCI allows for the quick switching between the high resolution continuous connectivity to Atlas-level discrete connectivity using the `parcellate_fc` and `parcellate_sc` functions. This means that it is possible to perform conventional analysis with any number of atlases, without the need to rerprocess data.

```
load('./example_data/example_fc.mat')
load('./example_data/example_sc.mat')

% load required SBCI data for mapping and analysis
[sbci_parc, sbci_mapping, adjacency] = load_sbci_data('./example_data/', 0.94);
 
% generate discrete connectivity 
dct_fc = parcellate_fc(fc, sbci_parc(3), sbci_mapping, 'roi_mask', [1,36]);
dct_sc = parcellate_sc(sc, sbci_parc(3), sbci_mapping, 'roi_mask', [1,36]);
```

![discretefig](https://user-images.githubusercontent.com/49790825/127187831-2e3ee7d9-4efc-4f75-af5f-031971b2ca4d.png)

For a complete example, see `example_discrete_connectivity.m`
