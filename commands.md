---
title: Commands
layout: default
has_children: false
permalink: /commands/
---

#### Table of contents
1. [`bngal-build-nets`](#bngal-build-nets)
2. [`bngal-summarize-nets`](#bngal-summarize-nets)


# `bngal-build-nets`

***Build a network from input abundance data.***

```
bngal-build-nets --help
```

### **`--asv_table`**

* **(Required)** Path to taxonomic count table named by Silva138-style taxonomies. (i.e., `d__DOMAIN;p__PHYLUM;c__CLASS;o__ORDER;f__FAMILY;g__GENUS;s__SPECIES`). Ideally rarefied and filtered as necessary.
	* First column must be named `'sample-id'` and must contain unique identifiers.
	* Must be an absolute abundance table.
	* If the input table is collapsed higher than level 7 (ASV/OTU), be sure to specify the `--taxonomic_level` option accordingly.
	* A formatted example table can be found here ([download](https://journals.asm.org/doi/suppl/10.1128/aem.01682-23/suppl_file/aem.01682-23-s0006.csv))

### **`--metadata`**
* **(Required)** Sample metadata corresponding to asv_table. Must be a .CSV file with sample identifiers in a column named `sample-id.`
	* A formatted example table can be found here ([download](https://journals.asm.org/doi/suppl/10.1128/aem.01682-23/suppl_file/aem.01682-23-s0004.csv))


### **`--output`**
* *(Optional)* Output directory for networks and graphs.<br>*Default = `./bngal-results`*

### **`--correlation`**
* *(Optional)* Metric for pairwise comparisons. Can be one of `pearson` or `spearman`.<br>*Default = `spearman`*

### **`--transformation`**
* *(Optional)* Numeric transformation to apply to input data before correlation calculations. Can be one of `log10` or `NULL`.<br>*Default = `NULL`*

### **`--corr_columns`**
* *(Optional)* Metadata columns to include in pairwise correlation networks. Multiple columns may be provided given the following syntax: `'col1,col2'`.<br>*Default = `NULL`*

### **`--corr_cutoff`**
* *(Optional)* Absolute correlation coefficient (defined from `--correlation`) cutoff for pairwise comparisons.<br>*Default = `0.6`*

### **`--p_value`**
* *(Optional)* Maximum cutoff for p-values calculated from pairwise relationships.<br>*Default = `0.05`*

### **`--abun_cutoff`**
* *(Optional)* Relative abundance cutoff for taxa (values 0-1 accepted). Anything lower than this value is removed before network construction.<br>*Default = `0`*

### **`--cores`**
* *(Optional)* Number of CPUs to use. Can only parallelize on Mac or Linux OS.<br>*Default = `1`*

### **`--subnetworks`**
* *(Optional)* Metadata column by which to split data in order to create separate networks. If not provided, `bngal` will create a single network from the input ASV table.<br>*Default = `NULL`*

### **`--taxonomic_level`**
* *(Optional)* Taxonomic level at which to construct co-occurrence networks. Must be at the same level or above the input `--asv_table`. Can be one of `phylum`, `class`, `order`, `family`, `genus`, or `asv`<br>*Default = `asv`*

### **`--direction`**
* *(Optional)* Direction for `--abun-cutoff`. Can be one of `greaterThan` or `lessThan`.<br>*Default = `greaterThan`*

### **`--sign`**
* *(Optional)* Type of pairwise relationship to filter for network construction. Can be one of `positive` (`--corr_cutoff` > 0), `negative` (`--corr_cutoff` < 0), or `all` (`--corr_cutoff` > 0 and `--corr_cutoff` < 0).<br>*Default = `all`*

### **`--obs_threshold`**
* *(Optional)* 'Observational threshold'. Minimum number of unique observations required for a given pairwise relationship to be included in the network.<br>*Default = `5`*

### **`--graph_layout`**
* *(Optional)* Type of igraph layout for output network plots. Refer to the [igraph documentation](https://igraph.org/r/html/latest/layout_.html) for the full list of options.<br>*Default = `layout_nicely`*


# `bngal-summarize-nets`

***Summarize network statistics from `bngal-build-nets` output.***

```
bngal-summarize-nets --help
```

### **`--asv_table`**

* **(Required)** See [`asv_table`](#--asv_table)

### **`--metadata`**

* **(Required)** See [`metadata`](#--metadata)

### **`--output`**

* *(Optional)* See [`output`](#--output)

### **`--taxonomic_level`**
* *(Optional)* See [`taxonomic_level`](#--taxonomic_level)

### **`--subnetworks`**
* *(Optional)* See [`subnetworks`](#--subnetworks)

### **`--fill_ebc_by`**
* *(Optional)* Metadata column by which to fill EBC composition plots.<br>*Default = `NULL`*

### **`--interactive`**
* *(Optional)* Determines whether output EBC composition plots are exported as interactive HTMLs (`TRUE`) or static PDFs (`FALSE`).<br>*Default = `FALSE`*

### **`--cores`**
* *(Optional)* See [`cores`](#--cores)

### **`--query`**
* *(Optional)* A query string to construct co-occurrence plots for a specific taxon. Be sure to use the full taxonomic ID as appropriate for the given taxonomic level as noted in the `*taxa_spread.csv` output in the network-summaries subfolder. Multiple queries may be provided given space characters: `Archaea;Crenarchaeota Bacteria;Actinobacteriota`<br>*Default = `NULL`*

### **`--tip_shape_by`**
* *(Optional)* Define the shape of the summary barplot's dendrogram tips by a metadata column.<br>*Default = `NULL`*

### **`--skip_plotting`**
* *(Optional)* Skip the plotting of taxonomic barplots and EBC composition plots. Useful if you want to test multiple `--query` inputs on the same data.<br>*Default = `FALSE`*

