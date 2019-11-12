#prepare gmt for GSEA
options(stringsAsFactors = F)
library(tidyverse)
egg <- read.delim("eggnog_annotation.txt", sep = "\t")
#If you want KOs or GOs rather than eggNOG, change "egg_anno" to "GO" or "KEGG_module"
gterms <- egg %>% dplyr::select(GID = query_name, GENENAME = egg_anno) %>% na.omit()
gene2egg <- data.frame(GID = character(), egg = character())
for (row in 1:nrow(gterms)) {
  gene_terms <- str_split(gterms[row,"GENENAME"], ",", simplify = FALSE)[[1]]  
  gene_id <- gterms[row, "GID"][[1]]
  tmp <- data_frame(GID = rep(gene_id, length(gene_terms)),
                    egg = gene_terms)
  gene2egg <- rbind(gene2egg, tmp)
}
gene2egg_unique <- unique.data.frame(gene2egg)
#get eggNOG annotation from eggNOG website http://eggnog5.embl.de/
annotation <- read.delim('~/Downloads/e5.og_annotations.tsv', sep = "\t", header=F, quote="",check.names=F)
annotation <- annotation[-1]
annotation <- annotation[!duplicated(annotation$V2),]
write.table(annotation,"unique_annotation.txt", sep="\t", quote=F, col.names=F, row.names=F)#save for convinence
gmt <- data.frame(unlist(tapply(gene2egg_unique$GID,as.factor(gene2egg_unique$egg),function(x) paste(x,collapse ="\t"),simplify =F)))
#gmt$description <- rep("An egg annotation", nrow(gmt))
gmt$GID <- row.names(gmt)
gmt <- merge(gmt,annotation,by.x="GID",by.y="V2")
gmt$querys <- gmt$unlist.tapply.gene2egg_unique.GID..as.factor.gene2egg_unique.egg...
gmt <- gmt[,-c(2,3)]
#avoid duplicate
gmt <- gmt[!duplicated(gmt$GID),]
write.table(gmt, "egg.gmt", col.names = F, row.names = F, quote = F, sep = "\t")
