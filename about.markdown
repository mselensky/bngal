---
layout: page
title: About
permalink: /about/
---

**Welcome to the `bngal` wiki!**

*What is `bngal`?*

Biological Network Graph Analysis and Learning (`bngal`) is a package primarily written in R to create high-quality, complex correlation networks from microbial abundance data.

`bngal` can create correlation networks at each level of taxonomic classification (phylum to ASV) from a taxonomic count table to visualize complex co-occurrence substructures in the data via [edge betweenness clustering][ebc]. Numeric variables from a corresponding metadata table can be optionally included to explore environmental-taxonomic correlations. "Subcommunity networks" can also be created in parallel to explore different correlation patterns within a dataset in addition to a global comparison. For example, one may want to examine separate networks for the human skin, oral, and gut microbiomes from the same dataset, while also examining microbial co-occurrence patterns across the entire body. Another may want to do the same thing for subsurface environments that span distinct geological contexts. As such, microbial ecologists from a wide range of backgrounds may be interested in applying bngal to model microbial niche space in the habitats they study with network analysis!

[ebc]: https://igraph.org/r/doc/cluster_edge_betweenness.html