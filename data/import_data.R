# package internal bngal data

# phylum color scheme (based on phylogeny; tol = Tree of Life)
phylum_colors_tol <- read.csv("data/silva_key_curated.csv")
# default EBC coloring scheme
ebc_colors <- read.csv("data/ebc_colors.csv")
# alternate phylum color scheme (based on culture status)
phylum_colors <- read.csv("data/silva138_phylum_colors.csv")

usethis::use_data(phylum_colors_tol, phylum_colors, ebc_colors, internal=TRUE, overwrite = TRUE)
