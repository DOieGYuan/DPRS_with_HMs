library(circlize)
library(grid)
library(ComplexHeatmap)
library(vegan)
#read proteins expression data
bins <- read.delim("bin_avg_morethan10.txt", header = T, row.names = 1, sep = "\t")
#Heatmap color span
color <- colorRamp2(c(0, 0.5, 1, 1.5, 2), c("#104E8B","#4F94CD", "#EECFA1", "#CD2626", "#660000"))
#Sample tree
sample_dist <- vegdist(t(bins))
sample_tree <- hclust(sample_dist, method = "average")
#Species tree
bins_dist <- vegdist(bins)
bins_tree <- hclust(bins_dist, method = "average")
#bins_hellinger <- data.frame(decostand(t(bins), method = "hellinger"))
#bins_tree <- hclust(dist(t(bins_hellinger)), method = "average")
#row_annotation
tax <- read.delim("tax.txt", sep = "\t", header = F)
tax = tax[, -2]
#if no ordered, do order step below first
ordered$Genomes <- row.names(ordered)
avg10_bins_anno <- merge(ordered, tax, by.x = "Genomes", by.y = "V1", all = F, sort = F)
write.table(avg10_bins_anno, "order_avgmorethan10_bins_annotation.txt", row.names = F, col.names = T, quote = F,sep = "\t")
#simplify the dataframe
avg10_anno <- read.delim("order_avgmorethan10_bins_annotation.txt", sep = "\t", row.names = 1)
row_anno <- HeatmapAnnotation(df = avg10_anno, 
                              col = list(Tax = c("Alphaproteobacteria"="#eee685",
                                                 "Betaproteobacteria"="#cdb5cd", 
                                                 "Gammaproteobacteria"="#ffb5c5", 
                                                 "Deltaproteobacteria"="#ff8c00", 
                                                 "Actinobacteria"="#7ccd7c",
                                                 "Flavobacteriia"="#c1ffc1", 
                                                 "Chitinophagia"="#cd4f39", 
                                                 "Cytophagia"="#4f94cd", 
                                                 "Nitrospira"="#aaaaaa", 
                                                 "Planctomycetes"="#7d26cd", 
                                                 "Verrucomicrobia"="#ff0066", 
                                                 "Gemmatimonadetes"="#008b00", 
                                                 "Bacteroidetes"="#33ccff")), 
                              na_col = "white", 
                              which = "row", 
                              show_legend = F)
#sample_annotation
col_anno <- data.frame(names(bins))
row.names(col_anno) <- col_anno$names.bins.
sample_anno <- HeatmapAnnotation(df = col_anno, 
                                 col = list(names.bins. = c("s1"="#c2c2c2",
                                                            "s2"="#c2c2c2", 
                                                            "s3"="#c2c2c2", 
                                                            "s4"="#0099ff", 
                                                            "s5"="#33cc66",
                                                            "s16"="#ffff66", 
                                                            "s18"="#0066cc", 
                                                            "s19"="#339933", 
                                                            "s20"="#ffcc33", 
                                                            "s22"="#003399", 
                                                            "s23"="#336600", 
                                                            "s24"="#ff9900")),
                                 show_annotation_name = F,
                                 which = "column",
                                 show_legend = F)
#plot heatmap
hm1 <- Heatmap(log10(bins), 
               col = color, 
               row_dend_width = unit(30, "mm"), 
               heatmap_legend_param = list(title= "lg(Abundance)", title_position ="topcenter", legend_direction="vertical", 
                                           row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 10)), 
               cluster_rows = bins_tree, 
               cluster_columns = sample_tree, 
               column_order = c("s24","s20","s16","s22","s23","s19","s5","s18","s4","s2","s3","s1"),
               column_dend_reorder = T, 
               show_row_names = F,
               top_annotation = sample_anno)
#get order
row_order <- row_order(hm1)
colum_order <- column_order(hm1)
ordered <- bins[row_order[[1]],colum_order]
#export
write.table(ordered,file="order_avgmorethan10_bins.txt",quote = FALSE,sep='\t')
#get functional bins
function_bins <- read.delim("fun_gene.txt", sep = "\t", row.names = 1)
row_func <- HeatmapAnnotation(df = ordered_func_gene,
                              na_col = "white", 
                              which = "row", 
                              name = "Functional genes", 
                              show_annotation_name = T, 
                              show_legend = F)
#plot density heatmap
densityHeatmap(log10(ordered))
#draw heatmap with row annotation
hm1 + row_anno + row_func
#
#
#plot fun_gene by ggplot2
ordered$Genome <- row.names(ordered)
function_bins$Genome <- row.names(function_bins)
ordered_func_gene <- merge(ordered, function_bins, by = "Genome", sort = F)[, -c(2:13)]
ordered_func_gene <- pivot_longer(ordered_func_gene,-Genome, names_to = "Gene", values_to = "PA")
ordered_func_gene[ordered_func_gene==""] <- NA
ordered_func_gene <- na.omit(ordered_func_gene)
ordered_func_gene$Gene <- factor(ordered_func_gene$Gene, levels = names(function_bins)[-12])
ordered_func_gene$Genome <- factor(ordered_func_gene$Genome, levels = ordered$Genome)
library(tidyverse)
ggplot(ordered_func_gene, aes(x = Gene, y = Genome)) + 
  geom_point(shape = 1, size = 1, color = "#473c8b") + 
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.border = element_blank()) + 
  xlab(NULL) +
  ylab(NULL)
#or heatmap
Heatmap(ordered_func_gene[-1],cluster_rows = F, cluster_columns = F, col = c("Yes"="#B22222", "NA"="#ffffff"))
