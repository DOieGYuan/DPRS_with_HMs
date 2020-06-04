setwd("/Performance")
library(tidyverse)
library(ggpubr)
# Significant analysis
pf <- read_tsv("TN_removal_rate.txt") %>% rename(time=Time)
cd <- filter(pf,metal == "Cd")
ni <- filter(pf,metal == "Ni")
cr <- filter(pf,metal == "Cr")
ck <- filter(pf,metal == "CK")
te <- rbind(cd,ck) %>% filter(time <= 48&time >=1)
te <- rbind(cd,ck) %>% filter(time <= 112&time >=56)
te <- rbind(cd,ck) %>% filter(time <= 186&time >=140)
te <- rbind(ni,ck) %>% filter(time <= 48&time >=1)
te <- rbind(ni,ck) %>% filter(time <= 106&time >=50)
te <- rbind(ni,ck) %>% filter(time <= 186&time >=122)
te <- rbind(cr,ck) %>% filter(time <= 48&time >=6)
te <- rbind(cr,ck) %>% filter(time <= 150&time >=98)
te <- rbind(cr,ck) %>% filter(time <= 242&time >=182)
te$metal <- factor(te$metal, levels = c("Cr", "CK"))
pdf("ttest_N.pdf", wi = 2, he = 2)
ggplot(te,aes(x = metal, y = tn)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1), 
                     expand = c(0,0)) +
  theme_bw() + 
  ylab("TN removal rate") + 
  stat_compare_means(method = "t.test",
                     paired = F,
                     label.x = 1,
                     label.y = 0.2,
                     method.args = list(alternative = "two.sided")) +
  theme(text = element_text(size = 4)) +
  xlab(NULL)
dev.off()
# P
pf <- read_tsv("P.txt")
cd <- filter(pf,metal == "Cd")
ni <- filter(pf,metal == "Ni")
cr <- filter(pf,metal == "Cr")
ck <- filter(pf,metal == "CK")
te <- rbind(cd,ck) %>% filter(time <= 48&time >=1)
te <- rbind(cd,ck) %>% filter(time <= 112&time >=56)
te <- rbind(cd,ck) %>% filter(time <= 186&time >=140)
te <- rbind(ni,ck) %>% filter(time <= 48&time >=1)
te <- rbind(ni,ck) %>% filter(time <= 106&time >=60)
te <- rbind(ni,ck) %>% filter(time <= 186&time >=122)
te <- rbind(cr,ck) %>% filter(time <= 48&time >=6)
te <- rbind(cr,ck) %>% filter(time <= 150&time >=98)
te <- rbind(cr,ck) %>% filter(time <= 236&time >=184)
te$metal <- factor(te$metal, levels = c("Cd", "CK")) # change as needed
pdf("ttest_P.pdf", wi = 2, he = 2)
ggplot(te,aes(x = metal, y = tp)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0,1), 
                     expand = c(0,0)) +
  theme_bw() + 
  ylab("TP removal rate") + 
  stat_compare_means(method = "t.test",
                     paired = T,
                     label.x = 1,
                     label.y = 0.2,
                     method.args = list(alternative = "two.sided")) +
  theme(text = element_text(size = 4)) +
  xlab(NULL)
dev.off()

# plot cols
dt <- read_tsv("Removal_rate_stablePhase.txt")
dt$stage <- factor(dt$stage,levels = c("LP","MP","HP"))
dt$metal <- factor(dt$metal,levels = c("Cd","Ni","Cr","CK_Cd","CK_Ni","CK_Cr"))
tn <- filter(dt,nutrient=="TN") %>% select(-nutrient)
tp <- filter(dt,nutrient=="TP") %>% select(-nutrient)

pdf("TP.pdf",wi=3.5,he=2)
ggplot(filter(tp,metal=="Cd"|metal=="Ni"|metal=="Cr"),
       aes(x=stage,y=tn)) +
  geom_col(fill="#d1d1d1") +
  geom_errorbar(aes(ymin=tn-se,ymax=tn+se), width=.3) +
  theme_pubclean() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     labels = scales::percent_format()) +
  facet_grid(.~metal,scales = "free_x") +
  theme(text = element_text(size=6))+
  guides(fill = F)
dev.off()

#combine
pdf("TNandTP.pdf",wi=4.5,he=2)
ggplot(filter(dt,metal=="Cd"|metal=="Ni"|metal=="Cr"),
       aes(x=stage,y=tn)) +
  geom_col(aes(fill=nutrient),position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=tn-se,ymax=tn+se), width=.15, 
                position = position_dodge(.9)) +
  theme_pubclean() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     labels = scales::percent_format()) +
  facet_grid(.~metal,scales = "free_x") +
  theme(text = element_text(size=6))+
  guides(fill = F) +
  scale_fill_brewer()
dev.off()
# line plot
dt$CK <- "Control"
dt$CK[dt$metal!="CK_Cd"&dt$metal!="CK_Ni"&dt$metal!="CK_Cr"]="Treated"
dt$metal[dt$metal=="CK_Cd"] <- "Cd"
dt$metal[dt$metal=="CK_Ni"] <- "Ni"
dt$metal[dt$metal=="CK_Cr"] <- "Cr"
pdf("6.TNandTP_line.pdf",wi=4.5,he=3)
ggplot(dt,
       aes(x=stage,y=tn,shape=CK,group=group)) +
  geom_point(aes(color=nutrient),size=2) +
  geom_line(aes(color=nutrient)) +
  geom_errorbar(aes(ymin=tn-se,ymax=tn+se,color=nutrient), width=.1) +
  theme_pubclean() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     labels = scales::percent_format()) +
  facet_grid(.~metal,scales = "free_x") +
  theme(text = element_text(size=6))+
  guides(fill = F) +
  scale_color_brewer()+
  scale_shape_manual(values = c("Control"=23,"Treated"=16))
dev.off()
# plot CK

pdf("TP_CK.pdf",wi=3.5,he=2)
ggplot(filter(tp,metal!="Cd"&metal!="Ni"&metal!="Cr"),
       aes(x=stage,y=tn)) +
  geom_point(shape=23) +
  theme_pubclean() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     labels = scales::percent_format()) +
  theme(text = element_text(size=6))+
  facet_grid(.~metal,scales = "free_x") +
  guides(fill = F)
dev.off()

pdf("TNandTP_CK.pdf",wi=4.5,he=2)
ggplot(filter(dt,metal!="Cd"&metal!="Ni"&metal!="Cr"),
       aes(x=stage,y=tn)) +
  geom_point(position = position_dodge(0.9),shape=18,size=2) +
  geom_errorbar(aes(ymin=tn-se,ymax=tn+se), width=.1) +
  theme_pubclean() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0),
                     labels = scales::percent_format()) +
  facet_grid(.~metal,scales = "free_x") +
  theme(text = element_text(size=6))+
  guides(fill = F) +
  scale_fill_brewer()
dev.off()
