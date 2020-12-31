### Plot Abundance profile of MAGs (Fig.2cd)
To create **Fig.2**, we use R package ggplot2.
```
setwd("/Taxonomy/Functional_microbes") # change to your working directory
library(tidyverse)
library(stringr)
info <- read_tsv("MAG_information.tsv")
data <- pivot_longer(info, -c(genome, AMO, HAO, LDA_socre, pvalue,
                              NXR, NAS, NAR, NAP, NIR, ccNIR, tax,
                              PPK1, PPK2,PPX),
                     names_to = "sample",
                     values_to = "abundance")
data$sample <- factor(data$sample, levels = c("Cd_O", "Cd_L", "Cd_M", "Cd_H",
                                              "Ni_O", "Ni_L", "Ni_M", "Ni_H",
                                              "Cr_O", "Cr_L", "Cr_M", "Cr_H", "CK"))
data <- mutate(data, name = str_c(genome, tax, sep = " "))
pick <- filter(data, NAR=="Yes"|NAP=="Yes"|NIR=="Yes") # plot potential denitrifiers (Fig. S2)
pick <- filter(data, NXR=="Yes") # plot potential nitrite oxidizers (Fig. S3)
pick <- filter(data, AMO=="Yes") # plot potential ammonia oxidizers Fig.2c (upper)
pick <- filter(data, AMO=="Yes"| genome == "MAG-46") # plot Fig.2c (lower)
pick <- filter(data, str_detect(ncbi.tax,"Candidatus Accumulibacter")|str_detect(ncbi.tax,"Dechloromonas"))
pick$abundance[pick$abundance==0] <- NA
pdf("figs/[please rename the picture name].pdf",
    wi = 3.5, he = 3)
    ggplot(pick, aes(x = sample, y = genome,
                 size = abundance, color = abundance)) +
  geom_point() +
  #scale_size(limits = c(0,400),breaks = c(0,100,200,400)) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.text.y.left = element_text(size = 3.5),
        legend.text = element_text(size = 4),
        axis.text.x.bottom = element_text(size = 4),
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  xlab(NULL) +
  ylab(NULL)
dev.off()
# functional profile of denitrifiers (Fig. S2 right panel)
func <- select(info, c(genome, NAP, NAR, NAS, NIR)) %>%
  filter(NAR=="Yes"|NAP=="Yes"|NIR=="Yes") %>%
  pivot_longer(-genome, values_to = "presence", names_to = "gene")
func$presence[func$presence=="Yes"] <- 1
pdf("figs/[renmae].pdf",
    wi = 1.5, he = 8)
ggplot(func, aes(x = gene, y = genome, color = presence)) +
  geom_point() +
  scale_size(limits = c(0,1),breaks = c(0, 1)) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.y.left = element_text(size = 3.5),
        legend.text = element_text(size = 4),
        axis.text.x.bottom = element_text(size = 4),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  xlab(NULL) +
  ylab(NULL)
dev.off()
# sum of functional taxa (Fig.2cd upper panel)
pdf("figs/[rename].pdf",
    wi = 9, he = 3)
ggplot(pick, aes(x = sample, y = abundance, fill = tax)) +
  geom_col() +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(limits = c(0,200),
                     expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3")
dev.off()
```
Note that we only show the ploting of microbes' abundance profiles.  
For ploting the performance of reactos in Fig.2ab, see [plot_performance.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/plot_performance.R) and [plot_Fig.S5.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/plot_Fig.S5.R).
