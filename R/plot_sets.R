plot_sets <- function(taxa.sets) {
  groups. <- taxa.sets %>%
    select(3:ncol(.)) %>%
    colnames(.)
  
  upset.plot <- ComplexUpset::upset(
    taxa.sets,
    groups.,
    base_annotations=list(
      'Intersection size'=intersection_size(
        counts=FALSE,
        mapping=aes(fill=phylum)) + 
        scale_fill_manual(values = phylum_color_dict) +
        theme(legend.position = "blank")
    ),
    width_ratio=0.1
  )
  
  upset.plot
}