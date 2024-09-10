---
layout: page
title: Quick start - global network analysis
permalink: /quick-start-global-net/
---

# Quick start guide: `bngal` 

# Global network analysis

Creating a "global network" means you are sending all of your taxonomic abundance data into the `bngal` pipeline. This is a useful starting analysis that visualizes broad trends across your entire dataset. The following global network example will analyze a [published dataset](https://journals.asm.org/doi/10.1128/aem.01682-23#supplementary-materials) to show the potential for regional biogeography of microbes living in the eastern Yucatan carbonate aquifer. After such a global analysis, we can [split up the data into separate networks](#separate-networks-by-metadata-column-region) based on a metadata column (`region`) to further refine pairwise correlation trends within each region.

## Step 1: `bngal-build-nets`

The first step in the bngal pipeline, `bngal-build-nets`, creates co-occurrence networks a specified level of taxonomic classification (phylum-ASV) and exports the output data for downstream processing. **Critically, the first column of the ASV/OTU table must be named `sample-id`, while the remaining columns are taxonomic IDs.** One of the metadata file's columns must also be named `sample-id` for subsequent mapping to the taxonomic abundance data (though column position does not matter in the metadata file). Both files must be in CSV format and contain unique `sample-id` values.

If you use [qiime2](https://qiime2.org/) to process your sequencing data like many microbial ecologists do, I recommend using the `read_qza()` function from the [qiime2R](https://github.com/jbisanz/qiime2R) package to import a [collapsed ASV-level table](https://docs.qiime2.org/2024.5/plugins/available/taxa/collapse/) into R and export it as a CSV file for use in bngal:

```
# import ASV table from qiime2
library(tidyverse)
library(qiime2R)
read_qza("collapsed-table-l7.qza") %>%
  .[["data"]] %>%
  t() %>%
  as.data.frame() %>%
  write_csv("example-asv-table.csv")
```

For the purposes of this tutorial, you can download Supplemental Tables S4 (a collapsed ASV-level table, rarefied to a depth of 9,957) and S2 (corresponding sample metadata) from [our paper](https://journals.asm.org/doi/10.1128/aem.01682-23#supplementary-materials) that first demonstrated `bngal`. In your terminal, you can run the following to create a folder named `bngal-test`, and go into it. 

Save the CSV files you downloaded into this folder, and name them variables called `META_DATA`, the sample metadata, and `TAX_TABLE`, the rarefied taxonomic abundance table. Finally, you can set an output directory, which we will call `all-communities` as we are going to start with creating one big network from the input dataset:

```
mkdir -p bngal-test && cd bngal-test

META_DATA=aem.01682-23-s0004.csv
TAX_TABLE=aem.01682-23-s0006.csv
OUT_DR=./all-communities
```

There are only three required options for `bngal-build-nets`: `--asv-table`, a rarefied ASV/OTU table, `--metadata`, sample metadata corresponding to `asv-table`, and `--output`, a directory path that must exist. By default, `bngal` will only create networks from pairwise associations that have at least 5 observations across the dataset and have an absolute correlation coefficient of at least 0.6 (`p` <= 0.05). `bngal` also assumes by defaults that co-occurrences will be analyzed at the ASV level from an ASV-level taxonomic count table, but users may tweak this and many other parameters to their liking - see the [`bngal-build-nets`](../commands/#bngal-build-nets) page for more details.

The simplest use case is to create a global network of the entire input ASV table without including any metadata variables. By default, the "observational threshold", or number of unique observations required per pairwise relationship to be included in the network, is set to 5. Building such a network looks like:

You might want to make such a "global" network when wanting to examine broad trends between groups of samples (communities) in your dataset. After we run `bngal-build-nets`, we'll move on to the next command, which will provide useful summaries and visualizations for us to interpret from the network that `bngal-build-nets` creates.

```
# using Docker
docker run -v `pwd`:/home/mambauser -it mjsel/bngal:dev \
  bngal-build-nets \
    --asv_table=$TAX_TABLE \
    --metadata=$META_DATA \
    --output=$OUT_DR
```

The above command results in several [output subfolders](../output_folders). The `pairwise-summaries` output subfolder contains a list of pairwise node statistics for each sample included in network analysis. The subfolder `network-plots` contains publication-ready network visualizations with nodes colored by a taxon's phylum and edge between cluster (EBC) in the `network-plots/pdfs` subfolder. In this example, nodes (individual ASV-level taxa) are colored by their EBC, a measure of how closely connected they are to other nodes in the network. The width of the edges (the lines connecting nodes) corresponds to a Spearman correlation co-efficient. Values can be positive (blue) or negative (red). As discussed in [our paper](https://journals.asm.org/doi/10.1128/aem.01682-23) (another shameless plug!), examining a heterogeneous dataset such as this one can result in a lot of smoothed relationships, but it's a good starting analysis. 

From this visualization, we can see that there are clear groups of microbes that tend co-occur more strongly than others; we are essentially seeing different cliques of 'crobes (EBCs) hanging out in the cafeteria (the entire dataset). One may (carefully) interpret some of these trends as the emergence of ecological niches; do some microbes co-occur with each other more because they share the same environmental/nutritional requirements? Do they outcompete others when nutrients are scarce, or is there even predation going on? Well, those are all great questions with pretty complicated answers, requiring a deeper look.

<object data="../figures/all-asv-by-edge_btwn_cluster.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="../figures/all-asv-by-edge_btwn_cluster.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="../figures/all-asv-by-edge_btwn_cluster.pdf">Download PDF</a>.</p>
    </embed>
</object>

Another way to visualize the network is to fill the node by microbial phylum (due to the size of the figure, the color key is saved to a separate PDF during runtime):

<object data="../figures/all-asv-by-phylum.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="../figures/all-asv-by-phylum.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="../figures/all-asv-by-phylum.pdf">Download PDF</a>.</p>
    </embed>
</object>

`bngal-build-nets` also provides these visualizations as interactive HTML plots, which is a great way to explore the networks in more detail! Hover over each node or edge for useful information. <a href="../figures/all-asv-network-coloredBy-edge_btwn_cluster.html" target="_blank">Click here to see the HTML</a> (it will open in a new tab).

Keep exploring the other `bngal-build-nets` [output files](../output_folders) as you'd like, but we are ready to move on to the next part of the pipeline.

## Step 2: `bngal-summarize-nets`

The second step in the bngal pipeline, `bngal-summarize-nets`, outputs more useful network summary data and plots. `bngal-summarize-nets` takes the output directory path of `bngal-build-nets` as its input. While `bngal-build-nets` constructs the networks and identifies edge betweenness clusters (EBCs) in the data, `bngal-summarize-nets` calculates the relative abundance of each EBC per sample in the dataset. These summary data, alongside the distribution of each EBC and taxon in the dataset, are exported to the `network-summary-tables` [subfolder](../output_folders). Notably, the `"*_tax_spread.csv"` output file reports the EBC assigned to a given taxon along with its abundance distribution in the data set.

`bngal-summarize-nets` is also useful to visualize biogeographic patterns of taxonomic and EBC distributions. For example, imagine that your samples are categorized by the metadata column `sample_type` and you want to examine whether certain EBCs are associated with certain types of samples. By including the `--fill_ebc_by` option below, `bngal-summarize-nets` will produce "EBC composition" plots that summarize which `sample_type` the majority of the taxa comprising each EBC originate. In this example dataset, several inferred hydrological "regions" were represented, so we wanted to see how the abundance of EBCs varied by each region to explore biogeography:

```
# using Docker
docker run -v `pwd`:/home/mambauser -it mjsel/bngal:dev \
  bngal-summarize-nets \
    --asv_table=$TAX_TABLE \
    --metadata=$META_DATA \
    --network_dir=$OUT_DR \
    --fill_ebc_by="region"
```

### Biogeographic trends of corresponding taxa

Running `bngal-summarize-nets` will produce a few more interesting figures to summarize the network. Since we specified `--fill_ebc_by="region"`, we can examine the overall composition of each EBC as it relates to the metadata column `region`. In other words, for a given EBC, such a composition reflects the total relative abundance of its taxa across different hydrological regions that were sampled. If taxa were evenly distributed across samples regardless of `region`, we might expect the filled values in the barplot to be the same. However, right away we see some clear trends that this is not the case; for example, taxa from EBC 10 is mostly found in `region4`, whereas those from EBC 20 are primarily mapped to `region6`:

<object data="../figures/asv-filled.by-region.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="../figures/asv-filled.by-region.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="../figures/asv-filled.by-region.pdf">Download PDF</a>.</p>
    </embed>
</object>

Though the above plot might reflect the composition of the edge betweenness clusters (EBCs) themselves, you are likely more interested in the composition of your samples with how they relate to EBCs. As discussed in [our paper](https://journals.asm.org/doi/10.1128/aem.01682-23#supplementary-materials) (yikes, how many plugs is that now?), an "EBC" can be thought of a microbial niche. The taxa within each EBC are shown to have [significant](../commands/#--p_value) correlative relationships with each other; those might be positive relationships (i.e., the relative abundance of one taxa increases with another), or they might be negative (one decreases as the other increases). `bngal` facilitates your data analysis by examining the abundance of each EBC across your dataset, while also clustering communities together based on their taxonomic composition:


<object data="../figures/asv-all-clustered-barplot-ebc.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="../figures/asv-all-clustered-barplot-ebc.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="../figures/asv-all-clustered-barplot-ebc.pdf">Download PDF</a>.</p>
    </embed>
</object>

Because we chose `--fill_ebc_by="region"`, the top hierarchical cluster is filled by the metadata column `region`. Generally speaking, in this dataset, microbial communities from a given hydrological region tend to be more similar to each other. We can clearly see this without even considering the EBC barplot below the hierarchical cluster. That said, when we begin to layer in the compositional data as it relates to EBCs, we see the emergence of potential niches that are perhaps defined by the `region`. For example, EBC 10 is most abundant in `region4`, while EBC 20 is most abundant in `EBC 20`. On further analysis, it turns out that the taxa mapped to these EBCs might be expected to have similar growth requirements that are consistent with the environment from which they are sampled. We also observe more ubiquitous EBCs, like EBC 3, that seem to be abundant in most regions. Biogeography is cool!

You can take a deep look into which individual taxa comprise a given EBC by examining the `network-summary-tables` output folder. Please see [this table](../output_folders/#output-subfolders) for a description of each file that is saved there. For more granular analysis of potential community structure, these files are very helpful.

As with every compositional and correlative analysis, it is up to you to decide if the relationships are actually meaningful *in situ*. But hopefully you find `bngal` to be a useful tool to help guide your analysis!

Now you might be wondering: if hydrological region is so important, might I observe deeper trends if I analyze separate networks for each region? The answer is yes, and `bngal` provides an easy way to do such an analysis. This is where the real fun begins. For a similar tutorial for multi-network analysis with `bngal`, please go on to the [next page](../multi-network_analysis).
