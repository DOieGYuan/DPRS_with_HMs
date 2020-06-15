setwd("D:/OneDrive/Study/paper/Revealing the resisitence mechanism of DPRS system/Taxonomy/Binning/CoNet")
library(tidyverse)
library(ggpubr)
fdr <- "FDR0.05"
nodes <- read_csv(paste("5",fdr,"node.csv",sep = "."))
edges <- read_csv(paste("5",fdr,"edge.csv",sep = "."))
pick <- filter(nodes, degree >= 10) %>%
  arrange(-degree, negdegree)
pick$genome <- factor(pick$genome, levels = pick$genome)
pick <-  select(pick, c(genome, posdegree, negdegree)) %>%
  pivot_longer(-genome, values_to = "degree", 
               names_to = "pn")
# plot
pdf(paste("6",fdr,"pdf",sep = "."), wi = 6, he = 3)
ggplot(pick, aes(x = genome, y=degree, fill =pn)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 4, 
                                          angle = 60, 
                                          vjust = 1,
                                          hjust =1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab(NULL)
dev.off()
# degree topology
info <- read_tsv("1.MAGinfo.txt")
info <- read_tsv("2.info_MAG.txt")#new
info1 <- left_join(nodes,info)
info1$func_tax[info1$func_tax!="PAO"] <- "Others"
ggplot(info1, aes(x = func_tax, y = posdegree, color = func_tax)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = .4) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  stat_compare_means()
# calculate topological properties
library(igraph)
node <- select(nodes, genome)
edge <- filter(edges,interactionType=="mutualExclusion") %>% select(c(Source, Target))
net <- graph_from_data_frame(d = edge, 
                             vertices = node, 
                             directed = F)
topology <- data.frame(degree(net, mode = "all", loops = F, normalized = F))
names(topology)[1] <- "degree"
edge_density(net)
topology$transitivity <- transitivity(net, type="local")
topology$betweeness <- betweenness(net, 
                                   directed = F, 
                                   weights = NULL, 
                                   nobigint = T,
                                   normalized = F)
topology$closeness <- closeness(net, mode = "all", 
                                weights = NULL,
                                normalized = F)
write.table(topology, paste("6",fdr,"topology.tsv",sep="."), 
            quote = F, sep = "\t")
topology <- read_tsv("6.FDR0.05.topology.tsv")
topology <- left_join(topology,info1)
ggplot(topology, aes(x = func_tax, y = betweeness, color = func_tax)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), alpha = .4) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1)) +
  stat_compare_means()

topology <- read_tsv("6.FDR0.05.topology.tsv")
topology.copresence <- read_tsv("6.FDR0.05.topology.copresence.tsv")
topology.mutualExclusion <- read_tsv("6.FDR0.05.topology.mutualExclusion.tsv")
top <- left_join(topology,topology.copresence)
top <- left_join(top,topology.mutualExclusion)
# plot degree
deg <- select(top, c(genome, degree,
                          posdegree, negdegree)) %>%
  pivot_longer(-genome, names_to = "pn", values_to = "degree")
deg$pn <- factor(deg$pn, levels = c("degree", 
                                    "posdegree",
                                    "negdegree"))
pdf(paste("6",fdr,"degree.distribution.pdf",sep = "."),
    wi = 3, he = 2.8)
ggplot(deg, aes(x = pn, y = degree, color = pn)) + 
  geom_violin(width = 0.6, scale = "width", trim = F) +
  #geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(width=0.1,position=position_dodge(0.9)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_color_brewer(palette = "Set2") +
  xlab(NULL) + 
  ylab("Degree") +
  scale_y_continuous(limits = c(0,30)) +
  guides(color=F)
dev.off()
# plot betw
bew <- select(top, c(genome, betweeness,
                     co.betweeness, ex.betweeness)) %>%
  pivot_longer(-genome, names_to = "pn", values_to = "bew")
pdf(paste("6",fdr,"betweenness.distribution.pdf",sep = "."),
    wi = 3, he = 2.8)
ggplot(bew, aes(x = pn, y = bew, color = pn)) + 
  geom_violin(width = 0.6, scale = "width", trim = F) +
  #geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(width=0.1,position=position_dodge(0.9)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_color_brewer(palette = "Set2") +
  xlab(NULL) + 
  ylab("Betweenness") +
  guides(color=F) +
  scale_y_log10()
dev.off()
# plot degree~betweenness
top <- left_join(top,info) %>%
  mutate(mean_ab = (s1+s2+s3+s4+s5+s16+s18+s19+s20+s22+s23+s24+s25)/13)
top$phylum[top$phylum=="Candidatus Hydrogenedentes"|
             top$phylum=="Cyanobacteria"|
             top$phylum=="Chloroflexi"|
             top$phylum=="Gemmatimonadetes"] <- "Others"
pdf(paste("6",fdr,"keystone.size.pdf",sep = "."),
    wi = 5, he = 2.5)
ggplot(top, aes(x = degree, y = betweeness, 
                color = phylum, size = mean_ab)) + 
  geom_point() +
  scale_size_continuous(range = c(.05, 4)) +
  theme_bw() +
  theme(legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),) +
  scale_color_brewer(palette = "Set3") +
  scale_y_sqrt()
dev.off()
# plot degree rate
keystone <- filter(top, degree >= 12, betweeness <= 200)
dr <- keystone %>%
  mutate(posrate = posdegree/degree) %>%
  select(genome, posrate) %>%
  arrange(-posrate)
dr$genome <- factor(dr$genome, levels = dr$genome)
dr <- dr %>%  mutate(negrate= 1- posrate) %>%
  pivot_longer(-genome, names_to = "pn", 
               values_to = "degreerate")
pdf(paste("6",fdr,"pos.rate.pdf",sep = "."),
    wi = 6, he = 2.8)
ggplot(dr, aes(x = genome, y = degreerate, fill = pn)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 3.5, 
                                          angle = 60, 
                                          vjust = 1,
                                          hjust =1))+
  guides(fill = F)
dev.off()
