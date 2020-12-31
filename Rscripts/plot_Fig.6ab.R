library(tidyverse)
library(ggpubr)
library(ggsci)

nodes <- read_tsv("node.tsv") # exported from Cytoscape
edges <- read_tsv("edge.tsv") # exported from Cytoscape

# calculate topological properties
library(igraph)
node = select(nodes, Node)
# Full network
edge = edges %>%
  select(c(Source, Target))
net <- graph_from_data_frame(d = edge,
                             vertices = node,
                             directed = F)
topology = data.frame(degree(net, mode = "all", loops = F, normalized = F))
topology$genome = row.names(topology)
names(topology)[1] <- "degree"
topology = as_tibble(topology) %>%
  select(genome,degree)
topology$transitivity = transitivity(net, type="local")
topology$betweeness = betweenness(net,
                                   directed = F,
                                   weights = NULL,
                                   nobigint = T,
                                   normalized = F)
topology$closeness = closeness(net, mode = "all",
                                weights = NULL,
                                normalized = F)
write_tsv(topology, "1.topology.tsv")
# Positive network
edge = edges %>%
  filter(interactionType=="copresence") %>%
  select(c(Source, Target))
net <- graph_from_data_frame(d = edge,
                             vertices = node,
                             directed = F)
topology = data.frame(degree(net, mode = "all", loops = F, normalized = F))
topology$genome = row.names(topology)
names(topology)[1] <- "degree"
topology = as_tibble(topology) %>%
  select(genome,degree)
topology$transitivity = transitivity(net, type="local")
topology$betweeness = betweenness(net,
                                  directed = F,
                                  weights = NULL,
                                  nobigint = T,
                                  normalized = F)
topology$closeness = closeness(net, mode = "all",
                               weights = NULL,
                               normalized = F)
write_tsv(topology %>%
            rename(posdegree=degree,
                   co.transitivity=transitivity,
                   co.betweeness=betweeness,
                   co.closeness=closeness),
          "1.topology.Copresence.tsv")
# Negative network
edge = edges %>%
  filter(interactionType=="mutualExclusion") %>%
  select(c(Source, Target))
net <- graph_from_data_frame(d = edge,
                             vertices = node,
                             directed = F)
topology = data.frame(degree(net, mode = "all", loops = F, normalized = F))
topology$genome = row.names(topology)
names(topology)[1] <- "degree"
topology = as_tibble(topology) %>%
  select(genome,degree)
topology$transitivity = transitivity(net, type="local")
topology$betweeness = betweenness(net,
                                  directed = F,
                                  weights = NULL,
                                  nobigint = T,
                                  normalized = F)
topology$closeness = closeness(net, mode = "all",
                               weights = NULL,
                               normalized = F)
write_tsv(topology %>%
            rename(negdegree=degree,
                   ex.transitivity=transitivity,
                   ex.betweeness=betweeness,
                   ex.closeness=closeness),
          "1.topology.MutualExclusion.tsv")

topology = read_tsv("1.topology.tsv")
topology.copresence <- read_tsv("1.topology.Copresence.tsv")
topology.mutualExclusion <- read_tsv("1.topology.MutualExclusion.tsv")
top <- left_join(topology,topology.copresence) %>%
  left_join(topology.mutualExclusion)
# plot degree
deg <- select(top, c(genome, degree,
                          posdegree, negdegree)) %>%
  pivot_longer(-genome, names_to = "pn", values_to = "degree")
deg$pn <- factor(deg$pn, levels = c("degree",
                                    "posdegree",
                                    "negdegree"))

ggplot(deg, aes(x = pn, y = degree, color = pn)) +
  geom_violin(width = 0.6, scale = "width", trim = F) +
  #geom_point(position = position_jitter(width = 0.2)) +
  geom_boxplot(width=0.1,position=position_dodge(0.9)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_color_nejm() +
  xlab(NULL) +
  ylab("Degree") +
  scale_y_continuous(limits = c(0,50),expand = c(0.02,0.02)) +
  guides(color=F)
ggsave("Fig.6a",wi = 3, he = 2.8)

# plot degree~betweenness
top = left_join(top,read_tsv("MAG_information_summary")) %>%
  mutate(mean_ab = (Cd_O+Ni_O+Cr_O+Cd_L_1+Cd_L_2+Cd_L_3+
                    Cd_M_1+Cd_M_2+Cd_M_3+Cd_H_1+Cd_H_2+
                    Cd_H_3+Ni_L_1+Ni_L_2+Ni_L_3+
                    Ni_M_1+Ni_M_2+Ni_M_3+Ni_H_1+Ni_H_2+
                    Ni_H_3+Cr_L_1+Cr_L_2+Cr_L_3+
                    Cr_M_1+Cr_M_2+Cr_M_3+Cr_H_1+Cr_H_2+
                    Cr_H_3+CK_1+CK_2+CK_3)/33,
         tax = case_when(str_detect(ncbi.tax,"Accumulibacter")~"Ca.Accumulibacter",
                         str_detect(ncbi.tax,"Dechloromonas")~"Dechloromonas",
                         str_detect(ncbi.tax,"Pseudoxanthomonas")~"Pseudoxanthomonas",
                         str_detect(ncbi.tax,"Nitrospira")~"Nitrospira",
                         str_detect(ncbi.tax,"Phreatobacter")~"Phreatobacter",
                         str_detect(ncbi.tax,"g__Zoogloea")~"Zoogloea",
                         str_detect(ncbi.tax,"Methyloversatilis")~"Methyloversatilis",
                         str_detect(ncbi.tax,"Sphingopyxis")~"Sphingopyxis",
                         str_detect(ncbi.tax,"Rubrivivax")~"Rubrivivax",
                         str_detect(ncbi.tax,"Comamonas")~"Comamonas",
                         str_detect(ncbi.tax,"Nitrosomonas")~"Nitrosomonas",
                         str_detect(ncbi.tax,"Competibacter")~"Competibacter",
                         str_detect(ncbi.tax,"Thauera")~"Thauera",
                         T~"Others"))

ggplot(top, aes(x = degree, y = betweeness,
                color = tax, size = mean_ab)) +
  geom_point() +
  scale_size_continuous(range = c(.2, 6)) +
  theme_bw() +
  theme(legend.key.height = unit(5, "pt"),
        legend.key.width = unit(3, "pt"),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  scale_color_d3(palette = "category20") +
  scale_y_sqrt()
ggsave("Fig.6b",wi = 5, he = 3)

keystone <- filter(top, degree >= 15, betweeness <= 250)
dr <- keystone %>%
  mutate(posrate = posdegree/degree) %>%
  select(genome, posrate) %>%
  arrange(-posrate)
dr$genome <- factor(dr$genome, levels = dr$genome)
dr <- dr %>%  mutate(negrate= 1- posrate) %>%
  pivot_longer(-genome, names_to = "pn",
               values_to = "degreerate")

ggplot(dr, aes(x = genome, y = degreerate, fill = pn)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(size = 3.5,
                                          angle = 60,
                                          vjust = 1,
                                          hjust =1))+
  guides(fill = F) +
  scale_fill_nejm()
ggsave("Keystone.posdegree.rate.pdf", wi = 6, he = 2.8)
