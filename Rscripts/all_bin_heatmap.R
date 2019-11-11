library(circlize)
library(grid)
library(ComplexHeatmap)
library(vegan)
#read proteins expression data
bins <- read.delim("all_bin.txt", header = T, row.names = 1, sep = "\t")
#Heatmap color span
color_absolute <- colorRamp2(c(0, 50, 100, 150, 200), c("#104E8B","#4F94CD", "#EECFA1", "#CD2626", "#660000"))
color <- colorRamp2(c(0, 0.5, 1, 1.5, 2), c("#104E8B","#4F94CD", "#EECFA1", "#CD2626", "#660000"))
#color <- colorRamp2(c(0, 0.5, 1, 1.5, 2), c("#6959cd","#8470ff", "#EECFA1", "#ee6363", "#ee4000"))
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
ordered$Genomes <- row.names(ordered)
all_bins_anno <- merge(ordered, tax, by.x = "Genomes", by.y = "V1", all = F, sort = F)
write.table(all_bins_anno, "all_bins_annotation.txt", row.names = F, col.names = T, quote = F,sep = "\t")
#simplify the dataframe
all_bins_anno <- read.delim("all_bins_annotation.txt", sep = "\t", row.names = 1)
row_anno <- HeatmapAnnotation(df = all_bins_anno, 
                              col = list(Tax = c("Alphaproteobacteria"="#eee685",
                                                 "Betaproteobacteria"="#cdb5cd", 
                                                 "Gammaproteobacteria"="#ffb5c5", 
                                                 "Deltaproteobacteria"="#ff8c00", 
                                                 "Actinobacteria"="#7ccd7c",
                                                 "Flavobacteriia"="#c1ffc1", 
                                                 "Chitinophagia"="#cd4f39", 
                                                 "Cytophagia"= "#4f94cd", 
                                                 "Nitrospira"="#aaaaaa", 
                                                 "Planctomycetes"="#7d26cd", 
                                                 "Verrucomicrobia"="#ff0066", 
                                                 "Gemmatimonadetes"="#008b00", 
                                                 "Bacteroidetes"="#33ccff", 
                                                 "Firmicutes"="black",
                                                 "Acidobacteria"="#cc00ff",
                                                 "Oligoflexia"="#ff0000")), 
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
#plot density heatmap
densityHeatmap(log10(ordered[,-13]))
#draw heatmap with row annotation
hm1 + row_anno
