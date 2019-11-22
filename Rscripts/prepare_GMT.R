# prepare gmt for GSEA and g:profile
# use the interproscan versus eggmapper results
################
# interproscan #
################
library(tidyverse)
library(GO.db)
godb <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
# write_tsv(godb, "GOdb.txt")
go_ipr <- read_tsv("ipr_mp_GO_annotation.txt", col_names = F)
gene2go_ipr <- go_ipr %>% dplyr::select(Gene = X1, GOID = X14) %>% na.omit() %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,"\\|"))+1), seq = ""), sep = "\\|") %>% 
  gather(key = "X", value = "GOID", -Gene) %>% dplyr::select(Gene, GOID) %>% 
  na.omit() %>% base::unique()
# write_tsv(gene2go_ipr, "gene2go_ipr.txt")
go_annot_ipr <- gene2go_ipr %>% left_join(godb)
# write_tsv(go_annot, "go_annot.txt")
go2gene_ipr <- gene2go_ipr %>% group_by(GOID) %>% 
  summarise(Gene = str_c(Gene, collapse = "\t")) %>% 
  mutate(Count = str_count(Gene, "\t")+1) %>% arrange(desc(Count)) %>% 
  left_join(godb)
# write_tsv(go2gene_ipr, "go2gene_ipr.txt")
# get wego
wego_ipr <- gene2go_ipr %>% group_by(Gene) %>% 
  summarise(GOID = str_c(GOID, collapse = ",")) %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,","))+1), seq = ""), sep = ",")
write_tsv(wego_ipr, "wego_ipr.txt", col_names = FALSE, na = "")
#make gmt for EA
gmt_ipr <- dplyr::select(go2gene_ipr, GOID = GOID, Desc = TERM, Gene = Gene)
write.table(gmt_ipr, "GO_interproscan.gmt", col.names = F, quote = F, row.names = F, sep = "\t")

################
# eggNOGmapper #
################
go_emp <- read_tsv("emapper_mp_GO_annotation.txt", col_names = T)
# make gene2go
gene2go_emp <- go_emp %>% dplyr::select(Gene = query_name, GOID = GOs) %>% na.omit() %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,","))+1), seq = ""), sep = ",") %>% 
  gather(key = "X", value = "GOID", -Gene) %>% dplyr::select(Gene, GOID) %>% 
  na.omit() %>% base::unique()
go_annot_emp <- gene2go_emp %>% left_join(godb)
go2gene_emp <- gene2go_emp %>% group_by(GOID) %>% 
  summarise(Gene = str_c(Gene, collapse = "\t")) %>% 
  mutate(Count = str_count(Gene, "\t")+1) %>% arrange(desc(Count)) %>% 
  left_join(godb)
wego_emp <- gene2go_emp %>% group_by(Gene) %>% 
  summarise(GOID = str_c(GOID, collapse = ",")) %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,","))+1), seq = ""), sep = ",")
write_tsv(wego_emp, "wego_emp.txt", col_names = FALSE, na = "")
gmt_emp <- dplyr::select(go2gene_emp, GOID = GOID, Desc = TERM, Gene = Gene)
write.table(gmt_emp, "GO_emapper.gmt", col.names = F, quote = F, row.names = F, sep = "\t")

#####################
# combine above two #
#####################
gene2go_merged <- base::unique(rbind(gene2go_ipr, gene2go_emp))
wego_merged <- gene2go_merged %>% group_by(Gene) %>% 
  summarise(GOID = str_c(GOID, collapse = ",")) %>%
  separate(GOID, paste0("X", 1:(max(str_count(.$GOID,","))+1), seq = ""), sep = ",")
write_tsv(wego_merged, "wego_merged.txt", col_names = FALSE, na = "")
go2gene_merged <- gene2go_merged %>% group_by(GOID) %>% 
  summarise(Gene = str_c(Gene, collapse = "\t")) %>% 
  mutate(Count = str_count(Gene, "\t")+1) %>% arrange(desc(Count)) %>% 
  left_join(godb)
gmt_merged <- dplyr::select(go2gene_merged, GOID = GOID, Desc = TERM, Gene = Gene)
write.table(gmt_merged, "GO_merged.gmt", col.names = F, quote = F, row.names = F, sep = "\t")
#######################################
# another method (consistent results) #
#######################################
# slower than upper one
gterms <- go %>% dplyr::select(GID = X1, GO = X14) %>% na.omit()
gene2go_2 <- data.frame(GID = character(), GO = character())
for (row in 1:nrow(gterms)) {
  gene_terms <- str_split(gterms[row,"GO"], "\\|", simplify = FALSE)[[1]]  
  gene_id <- gterms[row, "GID"][[1]]
  tmp <- data_frame(GID = rep(gene_id, length(gene_terms)),
                    GO = gene_terms)
  gene2go_2 <- rbind(gene2go_2, tmp)
}
gene2GO_unique <- unique.data.frame(gene2go_2)
annotation <- read.delim('/media/linyuan/SSD1/Enrichment/unique_annotation.txt', sep = "\t", header=F, quote="",check.names=F)
gmt <- data.frame(unlist(tapply(gene2egg_unique$GID,as.factor(gene2egg_unique$egg),function(x) paste(x,collapse ="\t"),simplify =F)))
#gmt$description <- rep("An egg annotation", nrow(gmt))
gmt$GID <- row.names(gmt)
gmt <- unique(merge(gmt,annotation,by.x="GID",by.y="V1"))
gmt$querys <- gmt$unlist.tapply.gene2egg_unique.GID..as.factor.gene2egg_unique.egg...
gmt <- gmt[,-c(2,3)]
new_gmt <- gmt[!duplicated(gmt$GID),]
write.table(new_gmt, "wcgna_avg0.5.gmt", col.names = F, row.names = F, quote = F, sep = "\t")
