# DOieGYuan
library(tidyverse)
metal = "Cr" # selected from "Cd", "Ni" and "Cr"
cd.onto <- read_csv(paste("5", metal, "node.csv", sep = ".")) %>%
  select(Node, NES) %>%
  rename(Target = Node)
cd <- read_csv(paste("5", metal, "edge.csv", sep = ".")) %>%
  filter(Type == "Gene2MAG") %>%
  select(Source, Target) %>%
  left_join(cd.onto)
cd$NES[cd$NES>0] <- "Pos"
cd$NES[cd$NES<0] <- "Neg"
cd.genome <- group_by(cd, Source, NES) %>%
  summarise(n = n())
cd.rank <- data.frame(table(cd$Source)) %>% arrange(-Freq)
#names(cd.genome)[c(1,2)] <- c("genome", "Frequency")
cd.genome$Source <- factor(cd.genome$Source, 
                           levels = cd.rank$Var1)
pdf(paste(metal, "genomeCount.pdf", sep = "."), wi = 4, he = 2.5)
ggplot(cd.genome,aes(x = Source, y= n, fill = NES)) + 
  geom_col() + 
  theme_bw() +
  scale_fill_brewer(palette = "Dark2") +
  theme(text = element_text(size = 4),
        axis.text.x = element_text(angle=60,
                                   hjust = 1,
                                   vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Associations") +
  xlab(NULL)
dev.off()