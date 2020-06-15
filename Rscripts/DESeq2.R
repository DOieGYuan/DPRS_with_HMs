# DESeq2
setwd("your directory")

library(DESeq2)
vignette("DESeq2") #Highly recommonded to read

# load data [Separated counts matrix and group infromation]
cts <- as.matrix(read.table("metaproteome.abundance.txt", sep = "\t",header = T,
                 row.names = 1))
coldata <- read.table("group.txt", sep = "\t",header = T,
                      row.names = 1)

# cts and coldata should be consistent in sample order.
cts <- cts[, rownames(coldata)]
(all(rownames(coldata) %in% colnames(cts))) & (all(rownames(coldata) == colnames(cts)))
# should be TRUE

# Construct DESeq data matrix
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group)
dds

# Pre-filtering (see vignette for independent filtering)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Specifying the reference level
dds$group <- relevel(dds$group, ref = "CK")

# Differential expression analysis
dds <- DESeq(dds)
resCd <- results(dds, contrast=c("group","Cd","CK"), 
               alpha = 0.05) # FDR < 0.05
resNi <- results(dds, contrast=c("group","Ni","CK"), 
                 alpha = 0.05) # FDR < 0.05
resCr <- results(dds, contrast=c("group","Cr","CK"), 
                 alpha = 0.05) # FDR < 0.05
resCd
resNi
resCr
summary(resCd)
mcols(resCd)$description

# Shrinkage of effect size (LFC estimates)
resultsNames(dds)
resLFC_Cd <- lfcShrink(dds, coef="group_Cd_vs_CK", 
                    type="apeglm", parallel = T)
resLFC_Ni <- lfcShrink(dds, coef="group_Ni_vs_CK", 
                    type="apeglm", parallel = T)
resLFC_Cr <- lfcShrink(dds, coef="group_Cr_vs_CK", 
                    type="apeglm", parallel = T)
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1e2,1e6); ylim <- c(-10,10)
plotMA(resLFC_Cd, xlim=xlim, ylim=ylim, main="Cd_CK")
plotMA(resLFC_Ni, xlim=xlim, ylim=ylim, main="Ni_CK")
plotMA(resLFC_Cr, xlim=xlim, ylim=ylim, main="Cr_CK")

# Independent hypothesis weighting (IHW)
# for p value adjustment of DESeq2 results.
library("IHW")
resIHW_Cd <- results(dds, filterFun=ihw, 
                  contrast=c("group","Cd","CK"), 
                  alpha = 0.05)
resIHW_Ni <- results(dds, filterFun=ihw, 
                     contrast=c("group","Ni","CK"), 
                     alpha = 0.05)
resIHW_Cr <- results(dds, filterFun=ihw, 
                     contrast=c("group","Cr","CK"), 
                     alpha = 0.05)
summary(resIHW_Cd)
summary(resIHW_Ni)
summary(resIHW_Cr)
metadata(resIHW_Cd)$ihwResult
mcols(resIHW_Cd)$description

# plot volcano
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1e2,1e6); ylim <- c(-10,10)
plotMA(resIHW_Cd, xlim=xlim, ylim=ylim, main="Cd_CK")
plotMA(resIHW_Ni, xlim=xlim, ylim=ylim, main="Ni_CK")
plotMA(resIHW_Cr, xlim=xlim, ylim=ylim, main="Cr_CK")
# plot Counts
d <- plotCounts(dds, gene=which.min(resIHW_Cd$padj), intgroup="group", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Output
write.csv(as.data.frame(resCd), file="Cd_DESeq2.csv")
write.csv(as.data.frame(resNi), file="Ni_DESeq2.csv")
write.csv(as.data.frame(resCr), file="Cr_DESeq2.csv")
write.csv(as.data.frame(resIHW_Cd), file="Cd_DESeq2_IHW.csv")
write.csv(as.data.frame(resIHW_Ni), file="Ni_DESeq2_IHW.csv")
write.csv(as.data.frame(resIHW_Cr), file="Cr_DESeq2_IHW.csv")

# Count data transformations
# Extracting transformed values (vsd or rld)
vsd <- vst(dds, blind=FALSE) # Variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) # Regularized log transformation
head(assay(vsd), 3)
head(assay(rld), 3)

# Data quality assessment by sample clustering and visualization
# Heatmap of the count matrix
library("pheatmap")
# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# Principal component plot of the samples
pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()