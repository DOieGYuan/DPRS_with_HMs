library(ggvegan)
library(plyr)
species <- read.delim('species_selected.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE)
species <- data.frame(t(species))
distance <- vegdist(species, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(species) - 1), eig = TRUE)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = "text")
summary(pcoa)
pcoa$eig
point <- data.frame(pcoa$point)
species_coordinate <- wascores(pcoa$points[,1:2], species)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
sample_site$metal <- factor(sample_site$metal, levels = c('None', 'Cd','Ni','Cr'))
sample_site$stress <- factor(sample_site$stress, levels = c('Negligible', 'Low', 'Medium', 'High'))
group_border <- ddply(sample_site, 'metal', function(df) df[chull(df[[2]], df[[3]]), ])
#plot
ggplot(sample_site, aes(PCoA1, PCoA2)) + 
  theme_bw() +
  theme(panel.grid.major = element_line(color = '#D3D3D3', size = 0.1), 
        text = element_text(color = "black"), 
        panel.grid.minor = element_line(color = '#EBEBEB', size = 0.1)) + 
  #geom_vline(xintercept = 0, color = '#696969', size = 0.4) + 
  #geom_hline(yintercept = 0, color = '#696969', size = 0.4) + 
  #geom_polygon(data = group_border, aes(fill = metal)) + 
  geom_point(aes(color = metal, shape = stress), size = 2, alpha = 0.8) + 
  scale_shape_manual(values = c(15, 16, 17, 18)) +  
  scale_color_manual(values = c('#BDBDBD', '#4169E1', '#43CD80', '#DAA520')) + 
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) + 
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) #+ 
  #annotate('text', label = 'CK', x = 0.05, y = 0.15, size = 5, colour = '#C673FF') + 
  #annotate('text', label = 'Cd', x = -0.08, y = -0.125, size = 5, colour = '#73D5FF') + 
  #annotate('text', label = 'Ni', x = 0.28, y = -0.125, size = 5, colour = '#49C35A') + 
  #annotate('text', label = 'Cr', x = -0.18, y = 0, size = 5, colour = '#FF985C')
