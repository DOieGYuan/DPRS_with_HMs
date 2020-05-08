# reactor performance
library(tidyverse)
library(ggpubr)
#for TN removal rate
reactors_TN <- read_tsv("2.TN_removal_rate.txt")
reactors_TN$metal <- factor(reactors_TN$metal, levels = c("CK","Cd","Ni","Cr"))
reactors_TN_Cd <- filter(reactors_TN,metal == "Cd"|metal == "CK")
reactors_TN_Ni <- filter(reactors_TN,metal == "Ni"|metal == "CK")
reactors_TN_Cr <- filter(reactors_TN,metal == "Cr"|metal == "CK")
#plot
pdf("2.reactors_TN_Cd.pdf", width = 2.2, height = 2)
#pdf("2.reactors_TN_Ni.pdf", width = 2.2, height = 2)
#pdf("2.reactors_TN_Cr.pdf", width = 2.8, height = 2)
ggplot(reactors_TN_Cd, aes(x = Time, y = tn, shape = metal)) + 
  geom_point(size = .4) + 
  geom_vline(xintercept = c(49, 107), linetype = "dashed", color = "#4F4F4F") + 
  #geom_vline(xintercept = c(49, 153), linetype = "dashed", color = "#4F4F4F") + 
  scale_shape_manual(values = c(1, 19)) + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .2, .4, .6, .8, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 186), breaks = c(0, 32, 64, 96, 128, 160, 192, 224, 256), expand = c(0, 0)) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  theme(axis.text = element_text(color = "black")) + 
  ylab(NULL)
dev.off()

#for Nitrogen speciation
TN <- read_tsv("3.nitrogen.txt")
TN$metal <- factor(TN$metal, levels = c("Cd","Ni","Cr"))
TN_Cd <- filter(TN,metal == "Cd") %>% 
  pivot_longer(-c(time,metal),names_to = "type",values_to = "nitrogen")
TN_Ni <- filter(TN,metal == "Ni")%>% 
  pivot_longer(-c(time,metal),names_to = "type",values_to = "nitrogen")
TN_Cr <- filter(TN,metal == "Cr")%>% 
  pivot_longer(-c(time,metal),names_to = "type",values_to = "nitrogen")
TN_Cd$type <- factor(TN_Cd$type, levels = c("ammonia","nitrite","nitrate"))
TN_Ni$type <- factor(TN_Ni$type, levels = c("ammonia","nitrite","nitrate"))
TN_Cr$type <- factor(TN_Cr$type, levels = c("ammonia","nitrite","nitrate"))
#plot
#pdf("reactors_TP_CK.pdf", width = 2.73, height = 2)
pdf("3.nitrogen_Cd.pdf", width = 2.13, height = 2)
pdf("3.nitrogen_Ni.pdf", width = 2.13, height = 2)
pdf("3.nitrogen_Cr.pdf", width = 2.73, height = 2)
ggplot(TN_Cr, aes(x = time, y = nitrogen, fill = type)) + 
  geom_area() + 
  #geom_vline(xintercept = c(49, 107), linetype = "dashed", color = "#4F4F4F") + 
  geom_vline(xintercept = c(49, 153), linetype = "dashed", color = "#4F4F4F") + 
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 242),breaks = c(0, 32, 64, 96, 128, 160, 192, 224, 256), expand = c(0, 0)) + 
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  theme(axis.text = element_text(color = "black")) + 
  ylab(NULL)
dev.off()

#for P
TP <- read_tsv("4.P.txt")
TP$metal <- factor(TP$metal, levels = c("CK","Cd","Ni","Cr"))
TP_Cd <- filter(TP,metal == "Cd") %>%
  pivot_longer(-c(time,metal),names_to = "type",values_to = "phosphorus")
TP_Ni <- filter(TP,metal == "Ni")%>%
  pivot_longer(-c(time,metal),names_to = "type",values_to = "phosphorus")
TP_Cr <- filter(TP,metal == "Cr")%>%
  pivot_longer(-c(time,metal),names_to = "type",values_to = "phosphorus")
TP_Cd$type <- factor(TP_Cd$type, levels = c("tpin","tpout","po4"))
TP_Ni$type <- factor(TP_Ni$type, levels = c("tpin","tpout","po4"))
TP_Cr$type <- factor(TP_Cr$type, levels = c("tpin","tpout","po4"))
#plot
#pdf("reactors_PO4_CK.pdf", width = 2.73, height = 2)
pdf("4.P_Cd.pdf", width = 2.13, height = 2)
pdf("4.P_Ni.pdf", width = 2.13, height = 2)
pdf("4.P_Cr.pdf", width = 2.73, height = 2)
ggplot(TP_Cr, aes(x = time, y = phosphorus, color = type)) + 
  geom_point(size = .4) + 
  #geom_vline(xintercept = c(49, 107), linetype = "dashed", color = "#4F4F4F") + 
  geom_vline(xintercept = c(49, 153), linetype = "dashed", color = "#4F4F4F") + 
  scale_color_brewer(palette = "Set1") + 
  scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 242),breaks = c(0, 32, 64, 96, 128, 160, 192, 224, 256), expand = c(0, 0)) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "bottom", legend.title = element_blank()) + 
  theme(axis.text = element_text(color = "black")) + 
  ylab(NULL)
dev.off()