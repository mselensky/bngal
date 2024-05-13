---
title: Output folders
layout: default
has_children: false
permalink: /output_folders/
---

# `bngal-cli` output folders

After running the `bngal-build-nets` and `bngal-summarize-nets` executables, you will find several subfolders in your output directory.

### Output subfolders

| Subfolder | Parent function | Contents |
| :---: | :---: | :--- |
| **`network-data/`** | `bngal-build-nets` &nbsp; &nbsp; | Underlying quality-filtered network node and edge data for each level of taxonomic classification (with _.rds_ extension - for direct import into R) |
| **`pairwise-summaries/`** | `bngal-build-nets` | CSV files summarizing pairwise data for each sample included in network analysis. Columns contain the following information: <br/> * _input.nodes_ - Number of input nodes <br/> * _input.pairwise_ - Number of possible pairwise relationships from input data. <br/> * _post_obs_filt.nodes_ - Number of nodes passing the "observational threshold". <br/> * _post_obs_filt.pairwise_ - Number of pairwise relationships passing the "observational threshold". <br/> * _QC.nodes_ - Number of quality-controlled nodes included in network model after filtering based on correlation and _p_-value cutoffs. <br/> * _QC.pairwise_ - Number of pairwise relationships included in network model after filtering based on correlation and _p_-value cutoffs. |
| **`runtime-tables/`** | `bngal-build-nets` | Runtime information for network construction functions |
| **`network-plots/`** | `bngal-build-nets` | Interactive or static network visualizations for each level of taxonomic classification |
| **`ebc-composition-plots/`** | `bngal-summarize-nets` | Coreness centrality summary data for nodes in each edge betweenness cluster (EBC) plotted against EBC compositions. EBC compositions are calculated from original relative abundance data mapped to the metadata variable of interest defined by the `--fill_ebc_by` option. |
| **`network-summary-tables/`** | `bngal-summarize-nets` | Taxonomic and EBC summary data for each constructed network. Taxa that do not fall into a defined EBC are labeled `0` in the "edge_btwn_cluster" columns. Subdirectories are organized by level of taxonomic classification and contain the following for each network: <br/> * *_[community_name]_tax_spread.csv* - Max, min, mean (+ standard deviation), and median abundance values for each taxon in the network. <br/> * *_[community_name]_ebc_distribution.csv* - Max, min, mean (+ standard deviation), and median abundance values for each EBC in the network. <br/> * *_[community_name]_ebc_abundance_per_sample.csv* - A list of each EBC's relative abundance per sample. |
| **`taxa-barplots/`** | `bngal-summarize-nets` | Taxa barplots clustered at a given level of taxonomic classification. Subdirectories are organized by filled type: <br/> * _ebc/_ - Filled by edge betweenness cluster (EBC). <br/> * _grouping/_ - Filled by [family functional groupings](https://github.com/mselensky/bngal/blob/main/data/16S_families.csv) (for family level and below only). <br/> * _phylum/_ - Filled by taxonomic phylum. |
| **`connectivity_plots/`** | `bngal-summarize-nets` | Connection summary plot for a defined level of taxonomic classification, visualized at the phylum level. <br/> See table below for column descriptions of the associated .csv file output(s) in this folder. |

### Column descriptions of the `connectivity_plots/*connections_data.csv` output file for a given taxonomic level of classification.

| Column name | Description |
| :---: | :--- |
| `tot_xions` | The total number of direct connections (primary edges) associated with the given node. |
| `pos_xions` | The number of positive primary edges. |
| `inter_ebc_xions` | The number of primary edges connecting the node to a different EBC subcluster. |
| `inter_phy_xions` | The number of primary edges connecting the node to a different phylum. |
| `unique_inter_phy_xions` | The number of unique phyla connected by primary edges. |
| `unique_inter_ebc_xions` | The number of unique EBCs connected by primary edges. |
| `n_taxa_per_phy` | The number of taxa associated with the given node’s phylum in the network. |
| `total_xions_k2` | The summed number of primary edges from all nodes directly connected to the given node (secondary edges) |
| `inter_ebc_k2` | The number of secondary edges that connect to a different EBC subcluster. |
| `inter_phy_k2` | The number of secondary edges that connect to a different phylum. |