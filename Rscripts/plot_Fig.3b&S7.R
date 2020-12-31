setwd("/functional_genes")
library(tidyverse)
library(stringr)
# Load DESeq2 resulting FC files
Cd <- read_csv("Cd_DESeq2_IHW.csv") %>%
  select(c(X1,log2FoldChange,padj)) %>%
  rename(protein = X1, Cd_log2fc = log2FoldChange, Cd_padj = padj)
Cr <- read_csv("Cr_DESeq2_IHW.csv") %>%
  select(c(X1,log2FoldChange,padj)) %>%
  rename(protein = X1, Cr_log2fc = log2FoldChange, Cr_padj = padj)
Ni <- read_csv("Ni_DESeq2_IHW.csv") %>%
  select(c(X1,log2FoldChange,padj)) %>%
  rename(protein = X1, Ni_log2fc = log2FoldChange, Ni_padj = padj)
log2fc <- left_join(Cd,Ni) %>% left_join(Cr)
# Load KEGG category
brite <- read_tsv("ko00001.aligned.formatted.txt") %>%
  select(c(ko, L1, L2))
# Load genome mapping file
genome_file = "nitrifiers_genome.txt"
map <- read_tsv("aligned_proteins.tsv") %>%
  select(c(ID2,genome,eggNOG,egg_annotation,ipr2,ipr_annotation,ko,description)) %>%
  rename(ipr=ipr2,protein=ID2) %>%
  filter(AOB=="Yes") %>%
nit <- read_tsv(genome_file) %>% unique() %>%
  left_join(log2fc)
select <- is.na(nit$Cd_padj)&is.na(nit$Ni_padj)&is.na(nit$Cr_padj)|
  (abs(nit$Cd_log2fc)<1)&(abs(nit$Ni_log2fc)<1)&(abs(nit$Cr_log2fc)<1)|
  (nit$Cd_padj>=0.05)&(nit$Ni_padj>=0.05)&(nit$Cr_padj>=0.05)
nit <- nit[!select,]
nit$Cd_log2fc[nit$Cd_padj>0.05] <- NA
nit$Ni_log2fc[nit$Ni_padj>0.05] <- NA
nit$Cr_log2fc[nit$Cr_padj>0.05] <- NA
wid <- select(nit, -c(Cd_padj, Ni_padj, Cr_padj)) %>%
  pivot_wider(names_from = genome,
              values_from = c(Cd_log2fc, Ni_log2fc, Cr_log2fc))
wid <- left_join(wid, brite) %>% unique()
wid <- wid[!duplicated(wid$protein),]
write_tsv(wid, paste(genome_file, ".wider.txt", sep = ""))
# plot bubble
wid <- read_tsv(paste(genome_file, ".wider.combined.txt", sep = ""))
lon <- pivot_longer(wid, -c(No,eggNOG,egg_annotation,ipr,
                            ipr_annotation,ko,description,
                            L1,L2),
                    names_to = "group",
                    values_to = "foldchange") %>%
  mutate(metal = str_sub(group, start = 1, end = 2))
lon$group <- str_sub(lon$group, start = 11)
lon$metal <- factor(lon$metal, levels = c("Cd", "Ni", "Cr"))
lon$foldchange[abs(lon$foldchange)<1] <- NA
lon$foldchange[lon$foldchange > 8] <- 8
lon$foldchange[lon$foldchange < -8] <- (-8)
pdf(paste(genome_file, ".pdf", sep = ""), wi = 5, he = 4)
ggplot(lon, aes(x = group, y = as_factor(No),
                 size = abs(foldchange), color = foldchange)) +
  geom_point() +
  scale_colour_gradient2(low = "#104E8B",
                         high = "#660000") +
  theme_bw() +
  theme(axis.text.y.left = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text.x.bottom = element_text(size = 4),
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~metal)
dev.off()

#########
# PAOs  #
#########
genome_file = "PAOs_genomes.txt"
map <- read_tsv("genome_mapped_by_proteins.tsv") %>%
  select(c(ID2,genome,eggNOG,egg_annotation,ipr2,ipr_annotation,ko,description)) %>%
  rename(ipr=ipr2,protein=ID2) %>%
  filter(PAO=="Yes") %>%
  write_tsv(genome_file)
nit <- read_tsv(genome_file) %>% unique() %>%
  left_join(log2fc)
select <- is.na(nit$Cd_padj)&is.na(nit$Ni_padj)&is.na(nit$Cr_padj)|
  (abs(nit$Cd_log2fc)<1)&(abs(nit$Ni_log2fc)<1)&(abs(nit$Cr_log2fc)<1)|
  (nit$Cd_padj>=0.05)&(nit$Ni_padj>=0.05)&(nit$Cr_padj>=0.05)
nit <- nit[!select,]
nit$Cd_log2fc[nit$Cd_padj>0.05] <- NA
nit$Ni_log2fc[nit$Ni_padj>0.05] <- NA
nit$Cr_log2fc[nit$Cr_padj>0.05] <- NA
wid <- select(nit, -c(Cd_padj, Ni_padj, Cr_padj)) %>%
  pivot_wider(names_from = genome,
              values_from = c(Cd_log2fc, Ni_log2fc, Cr_log2fc))
wid <- left_join(wid, brite) %>% unique()
wid <- wid[!duplicated(wid$protein),]
write_tsv(wid, paste(genome_file, ".wider.txt", sep = ""))
# Manually combine identical functionality
wid <- read_tsv(paste(genome_file, ".wider.combined.txt", sep = "")) %>%
  filter(L2=="09101 Carbohydrate metabolism"|
           L2=="09102 Energy metabolism") %>%
  mutate(L2="09102 Energy metabolism") %>%
  left_join(read_tsv("ko00001.aligned.formatted.txt") %>%
              select(c(ko, L1, L2, L3)))
# write_tsv(wid, paste(genome_file, ".wider.txt", sep = ""))

# plot bubble
wid <- read_tsv(paste(genome_file, ".wider.txt", sep = ""))
lon <- pivot_longer(wid, -c(No,eggNOG,egg_annotation,ipr,
                            ipr_annotation,ko,description,
                            L1,L2,L3),
                    names_to = "group",
                    values_to = "foldchange") %>%
  mutate(metal = str_sub(group, start = 1, end = 2))
lon$group <- str_sub(lon$group, start = 11)
lon$metal <- factor(lon$metal, levels = c("Cd", "Ni", "Cr"))
lon$foldchange[abs(lon$foldchange) < 1] <- NA
lon$foldchange[lon$foldchange > 8] <- 8
lon$foldchange[lon$foldchange < -8] <- (-8)
pdf(paste(genome_file, ".pdf", sep = ""), wi = 7, he = 12)
ggplot(lon, aes(x = group, y = as_factor(No),
                size = abs(foldchange), color = foldchange)) +
  geom_point() +
  scale_colour_gradient2(low = "#104E8B",
                         high = "#660000") +
  theme_bw() +
  theme(axis.text.y.left = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.text.x.bottom = element_text(size = 4),
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~metal)
dev.off()
