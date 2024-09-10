---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
---

# Welcome to the `bngal` wiki!

What is `bngal`?

Biological Network Graph Analysis and Learning (`bngal`) is a package written in R to create high-quality, complex correlation networks from microbial abundance data.

<object data="https://mselensky.github.io/bngal/figures/all-asv-by-edge_btwn_cluster.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://mselensky.github.io/bngal/figures/all-asv-by-edge_btwn_cluster.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="https://mselensky.github.io/bngal/figures/all-asv-by-edge_btwn_cluster.pdf">Download PDF</a>.</p>
    </embed>
</object>

`bngal` can create correlation networks at each level of taxonomic classification (phylum to ASV) from a taxonomic count table to visualize complex co-occurrence substructures in the data via [edge betweenness clustering](https://igraph.org/r/doc/cluster_edge_betweenness.html). Numeric variables from a corresponding metadata table can be optionally included to explore environmental-taxonomic correlations. “Subcommunity networks” can also be created in parallel to explore different correlation patterns within a dataset in addition to a global comparison. For example, one may want to examine separate networks for the human skin, oral, and gut microbiomes from the same dataset, while also examining microbial co-occurrence patterns across the entire body. Another may want to do the same thing for subsurface environments that span distinct geological contexts. As such, microbial ecologists from a wide range of backgrounds may be interested in applying bngal to model microbial niche space in the habitats they study with network analysis!

Visit the [Quick Start Guide](quick-start-global-net) for examples on how to use `bngal`.

If you found `bngal` helpful for your research, please consider citing [our paper](https://doi.org/10.1128/aem.01682-23) that first used this tool.

Contact bngal.help@gmail.com with any questions, or [open a GitHub issue](https://github.com/mselensky/bngal/issues) if you have trouble with the software.

A huge thanks goes to Prof. [Maggie Osburn](https://www.earth.northwestern.edu/our-people/faculty/osburn-magdalena.html) and Dr. [Caitlin Casar](https://www.caitlincasar.com/) for their invaluable advice on the theory behind and implementation of `bngal`!

*`bngal` developer: Dr. [Matt Selensky](https://mselensky.github.io)*