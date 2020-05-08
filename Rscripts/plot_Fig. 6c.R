# DOieGYuan: plot Fig. 6c
library(tidyverse)
library(ggalluvial)
edge <- read_csv("5.FDR0.05.edge.csv") %>%
  select(Source,Target,interactionType)
tax <- read_tsv("2.info_MAG.txt") %>%
  select(genome,phylum,tax) %>%
  rename(Target=genome)
# find all interactions of specific genus and put them onto the first col
genus <- "Ca.Accumulibacter"
ls <- tibble(c("MAG-136", "MAG-163", "MAG-169", "MAG-105",
        "MAG-256", "MAG-20","MAG-28"))
genus <- "Dechloromonas"
ls <- tibble(c("MAG-137", "MAG-165", "MAG-189", "MAG-211",
               "MAG-44"))
genus <- "Pseudoxanthomonas"
ls <- tibble(c("MAG-154"))
genus <- "Zoogloea"
ls <- tibble(c("MAG-263"))
genus <- "Rubrivivax"
ls <- tibble(c("MAG-41","MAG-55","MAG-219","MAG-262"))
genus <- "Sphingopyxis"
ls <- tibble(c("MAG-104", "MAG-164"))
genus <- "Competibacter"
ls <- tibble(c("MAG-125", "MAG-192"))
genus <- "Nitrospira"
ls <- filter(tax, phylum=="Nitrospirae") %>% select(Target)
genus <- "Proteobacteria"
ls <- filter(tax, phylum==genus) %>% select(Target)
genus <- "Bacteroidetes"
ls <- filter(tax, phylum==genus) %>% select(Target)
genus <- "Planctomycetes"
ls <- filter(tax, phylum==genus) %>% select(Target)
genus <- "Taibaiella"
ls <- tibble(c("MAG-10", "MAG-2", "MAG-63", "MAG-72"))
genus <- "Planctomycetes_UBA2386"
ls <- tibble(c("MAG-90", "MAG-237", "MAG-184"))

names(ls) <- "Source"
ft1 <- left_join(ls,edge) %>% na.omit()
names(ls) <- "Target"
ft2 <- left_join(ls,edge) %>% select(Target, Source, interactionType) %>% na.omit()
names(ft2) <- c("Source", "Target", "interactionType")
ft <- rbind(ft1,ft2) %>% unique()
# taxonomy of "Target"
ft <- left_join(ft,tax)
# format for plot
dat <- ft %>% group_by(interactionType, phylum) %>%
  summarise(n=n())
# plot
pdf(paste("8",genus,"pdf",sep = "."),wi = 3,he = 2)
ggplot(dat, aes(y = n,
                axis1 = interactionType, axis2 = phylum)) +
  geom_alluvium(aes(fill = phylum), width = 1/12) +
  geom_stratum(width = 1/6) +
  #geom_label(stat = "stratum", label.strata = TRUE, label.size = .5) +
  scale_x_discrete(limits = c("Int", "Phylum"), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"))
dev.off()