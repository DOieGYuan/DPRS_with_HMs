# prepare gmt for GSEA and g:profile
# use the enrichM, interproscan and eggmapper results
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