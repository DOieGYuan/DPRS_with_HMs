#prepare gmt for GSEA
options(stringsAsFactors = F)
library(tidyverse)
egg <- read.delim("trimed.txt", sep = "\t")
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
gmt <- data.frame(unlist(tapply(gene2egg_unique$GID,as.factor(gene2egg_unique$egg),function(x) paste(x,collapse ="\t"),simplify =F)))
gmt$description <- rep("An egg annotation", nrow(gmt))
gmt$query_list <- gmt$unlist.tapply.gene2egg_unique.GID..as.factor.gene2egg_unique.egg...
gmt <- gmt[, -1]
write.table(gmt, "egg.gmt", col.names = F, row.names = T, quote = F, sep = "\t")
