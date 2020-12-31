# DOieGYuan: plot Fig. 6c
library(tidyverse)
library(ggalluvial)

edge = read_tsv("edge.tsv") %>%
  select(Source,Target,interactionType)
tax = read_tsv("MAG_information_summary.tsv") %>%
  select(genome,ncbi.tax) %>%
  mutate(phylum=str_remove(ncbi.tax,"d__Bacteria;p__")) %>%
  mutate(phylum=str_remove(phylum,";c__.*")) %>%
  mutate(genus=str_remove(ncbi.tax,".*f__[A-Z]*[a-z]* *[A-Z]*[a-z]*;g__")) %>%
  mutate(genus=str_remove(genus,";s__.*"))

# find all interactions of specific genus and put them onto the first col
interest.genus = c("Candidatus Accumulibacter",
                   "Dechloromonas",
                   "Pseudoxanthomonas",
                   "Methyloversatilis",
                   "Zoogloea",
                   "Rubrivivax",
                   "Sphingopyxis",
                   "Competibacter",
                   "Nitrospira",
                   "Taibaiella",
                   "Chthonomonas",
                   "Tabrizicola",
                   "Devosia",
                   "Candidatus Odyssella",
                   "Candidatus Obscuribacter")
for(now.genus in interest.genus){
  ls = filter(tax,genus==now.genus) %>%
    select(genome)

  names(ls) = "Source"
  ft1 = left_join(ls,edge) %>% na.omit()
  names(ls) = "Target"
  ft2 = left_join(ls,edge) %>% select(Target, Source, interactionType) %>% na.omit()
  names(ft2) <- c("Source", "Target", "interactionType")
  ft = rbind(ft1,ft2) %>% unique()
  # taxonomy of "Target"
  ft = left_join(ft,tax %>%
                   select(genome,phylum) %>%
                   rename(Target=genome))
  # format for plot
  dat = ft %>% group_by(interactionType, phylum) %>%
    summarise(n=n())
  # plot
  ggplot(dat, aes(y = n,
                  axis1 = interactionType, axis2 = phylum)) +
    geom_alluvium(aes(fill = phylum), width = 1/12) +
    geom_stratum(width = 1/6) +
    #geom_label(stat = "stratum", label.strata = TRUE, label.size = .5) +
    scale_x_discrete(limits = c("Int", "Phylum"), expand = c(0, 0)) +
    scale_fill_d3(palette = "category20c",drop=F) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(text = element_text(size = 6),
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(3, "pt"))
  ggsave("Fig.6c",device = "pdf",wi = 3,he = 2)
}
