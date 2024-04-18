# SBCI_Toolkit for SBCI Result Visualization

## Overview

The SBCI Toolkit is specifically designed for the visualization and analysis of structural and functional brain connectivity (SC/FC) data derived from the [SBCI Pipeline.](https://github.com/sbci-brain/SBCI_Pipeline) It offers robust tools to plot continuous SC/FC data, adjust resolution, compute Surface-Based Connectivity Integration (SBCI), and display various metrics on the cortical surface.

## Installation

This toolkit is implemented in MATLAB. Follow these steps for setup:

1. Clone or download the repository to your local machine.
2. Open MATLAB and navigate to the directory containing the toolkit.
3. Run the `main.m` script.

## Required Data

The toolkit is designed to analyze data outputs from the SBCI pipeline. Specifically, it requires two key files:

* `fc_avg_ico4.mat`: This file contains the functional connectivity (FC) matrix.
* `smoothed_sc_avg_0.005_ico4.mat`: This file contains the smoothed structural connectivity (SC) matrix.

Additionally, example data is available in the directory `example_data/SBCI_Individual_Subject_Outcome` for demonstration purposes.

## Key Parameters

### `atlas_index` Parameter

The `atlas_index` parameter specifies the index within the `sbci_parc` array that selects a particular brain parcellation for analysis. The `sbci_parc` array contains detailed information about each atlas, including its name and specific partitions. This setup allows users to access and utilize various brain parcellations supported by the toolkit.

Below is a list of available brain parcellations that can be accessed via the `atlas_index`:

| Index | Parcellation Name                        |
| ----- | ---------------------------------------- |
| 1     | BN_Atlas                                 |
| 2     | HCPMMP1                                  |
| 3     | PALS_B12_Brodmann                        |
| 4     | PALS_B12_Lobes                           |
| 5     | PALS_B12_OrbitoFrontal                   |
| 6     | PALS_B12_Visuotopic                      |
| 7     | Schaefer2018_1000Parcels_7Networks_order |
| 8     | Schaefer2018_100Parcels_7Networks_order  |
| 9     | Schaefer2018_200Parcels_7Networks_order  |
| 10    | Schaefer2018_300Parcels_7Networks_order  |
| 11    | Schaefer2018_400Parcels_7Networks_order  |
| 12    | Schaefer2018_500Parcels_7Networks_order  |
| 13    | Schaefer2018_600Parcels_7Networks_order  |
| 14    | Schaefer2018_700Parcels_7Networks_order  |
| 15    | Schaefer2018_800Parcels_7Networks_order  |
| 16    | Schaefer2018_900Parcels_7Networks_order  |
| 17    | Yeo2011_17Networks_N1000                 |
| 18    | Yeo2011_7Networks_N1000                  |
| 19    | aparc.a2005s                             |
| 20    | aparc.a2009s                             |
| 21    | aparc                                    |
| 22    | oasis.chubs                              |

To view and select a specific atlas for your analysis, you can adjust the `atlas_index` parameter accordingly. For example, to use the `aparc`:

```matlab
atlas_index = 21; % Example using 'aparc' atlas

disp(fprintf('the atlas current use is: %s', sbci_parc(atlas_index).atlas{1}));
```

### `roi_exclusion_index` Parameter

The `roi_exclusion_index` parameter allows users to selectively exclude specific regions from analysis and display. This feature is particularly useful for focusing on relevant areas and omitting regions that may not be of interest or are irrelevant to the current study. The `sbci_parc` array contains detailed information about each atlas, including its name and specific partitions.

To utilize this feature, define the `roi_exclusion_index` with the indices of the regions you wish to exclude from the matrix plot. Each index corresponds to a specific region in the atlas used. For example, if regions labeled as 'LH_missing' and 'RH_missing' are not required for your analysis, you can exclude them by setting:

```matlab
roi_exclusion_index = [1,36]; % Index 1: 'LH_missing', Index 36: 'RH_missing'
```

## Usage

The main functionality of the SBCI Toolkit is encapsulated in the `main.m` script. Example data are in the `example_data` folder. The `main.m` can do the following visualization and connectivity data analysis:

### **Display Continuous SC/FC Data**

Visualize high-resolution structural and functional connectivity (SC/FC) data. The script accepts structural and functional connectivity matrices, applies specified atlas and exclusion indices, and displays the connectivity patterns.

```matlab
plot_sbci_mat(fc, sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 1, 'clim', [-0.1, 0.35]);

plot_sbci_mat(log((10^7*sc) + 1), sbci_parc(atlas_index), 'roi_mask', roi_exclusion_index, 'figid', 2, 'clim', [0, 3.5]);

```

![](https://raw.githubusercontent.com/ytr1023/img/main/continuous_fcsc.png)

### **Convert Continuous SC/FC to Discrete ROI based SC/FC**

This functionality allows for the transformation of continuous connectivity data into discrete ROI-based SC/FC, which is suitable for more conventional atlas-based analysis. This process is handled without the need to reprocess the underlying data.

```matlab
adjust_connectivity_res_for_fc(fc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);

adjust_connectivity_res_for_sc(sc, sbci_parc, atlas_index, sbci_mapping, roi_exclusion_index);

```

![](https://raw.githubusercontent.com/ytr1023/img/main/continueToDiscrete.png)

### **Compute and Display SFC on the Cortical Surface**

Computes and visualizes both global and local surface functional connectivity (SFC) directly on the cortical surface. This function leverages connectivity data alongside surface mapping and parcellation details to provide a detailed view of connectivity patterns.

```matlab
sfc_gbl = calculate_sfc_gbl(sc, fc, 'triangular', false);
sfc_loc = calculate_sfc_loc(sc, fc, sbci_parc(atlas_index), 'triangular', false);

plot_cortical_sfc(sfc_gbl, sfc_loc, sbci_surf, sbci_mapping, sbci_parc, atlas_index);
```

![](https://raw.githubusercontent.com/ytr1023/img/main/sfc.png)

### **Display Scalar Value on the SBCI Cortical Surface**

Enables the visualization of scalar values such as genetic heritability or cortical thickness on the cortical surface. This feature utilizes the surface mapping and parcellation data to overlay scalar data effectively.

```matlab
plot_value_cortically(sbci_surf, sbci_mapping, 'example.txt');
```

![](https://raw.githubusercontent.com/ytr1023/img/main/value.png)
