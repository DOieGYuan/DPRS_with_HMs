setwd("D:/OneDrive/Study/paper/Revealing the resisitence mechanism of DPRS system/Revision1/MAG_taxonomy")

library(tidyverse)
# Belows derived from GTDBvsNCBI
phylum = read_tsv("0.gtdb2ncbi_phylum.txt") %>% 
  mutate(ncbi.phylum=case_when(ncbi.phylum=="p__"~str_c(gtdb.phylum,"*",sep=""),
                               T~ncbi.phylum))
write_tsv(phylum,"0.gtdb2ncbi_phylum_c.txt")

class = read_tsv("0.gtdb2ncbi_class.txt") %>% 
  mutate(ncbi.class=case_when(ncbi.class=="c__"~str_c(gtdb.class,"*",sep=""),
                               T~ncbi.class))
write_tsv(class,"0.gtdb2ncbi_class_c.txt")

order = read_tsv("0.gtdb2ncbi_order.txt") %>% 
  mutate(ncbi.order=case_when(ncbi.order=="o__"~str_c(gtdb.order,"*",sep=""),
                              T~ncbi.order))
write_tsv(order,"0.gtdb2ncbi_order_c.txt")

family = read_tsv("0.gtdb2ncbi_family.txt") %>% 
  mutate(ncbi.family=case_when(ncbi.family=="f__"~str_c(gtdb.family,"*",sep=""),
                            T~ncbi.family))
write_tsv(family,"0.gtdb2ncbi_family_c.txt")

genus = read_tsv("0.gtdb2ncbi_genus.txt") %>% 
  mutate(ncbi.genus=case_when(ncbi.genus=="g__"~str_c(gtdb.genus,"*",sep=""),
                            T~ncbi.genus))
write_tsv(genus,"0.gtdb2ncbi_genus_c.txt")

species = read_tsv("0.gtdb2ncbi_species.txt") %>% 
  mutate(ncbi.species=case_when(ncbi.species=="s__"~str_c(gtdb.species,"*",sep=""),
                            T~ncbi.species))
write_tsv(species,"0.gtdb2ncbi_species_c.txt")

xxx = read_tsv("0.gtdb2ncbi_xxx.txt") %>% 
  mutate(ncbi.xxx=case_when(ncbi.xxx=="c__"~str_c(gtdb.xxx,"*",sep=""),
                            T~ncbi.xxx))
write_tsv(xxx,"0.gtdb2ncbi_xxx_c.txt")

# Belows are NCBIvsGTDB
phylum = read_tsv("0.ncbi2gtdb_phylum.txt") %>% 
  mutate(gtdb.phylum=case_when(gtdb.phylum=="p__"~str_c(ncbi.phylum,"*",sep=""),
                               T~gtdb.phylum))
write_tsv(phylum,"0.ncbi2gtdb_phylum_c.txt")

class = read_tsv("0.ncbi2gtdb_class.txt") %>% 
  mutate(gtdb.class=case_when(gtdb.class=="c__"~str_c(ncbi.class,"*",sep=""),
                              T~gtdb.class))
write_tsv(class,"0.ncbi2gtdb_class_c.txt")

order = read_tsv("0.ncbi2gtdb_order.txt") %>% 
  mutate(gtdb.order=case_when(gtdb.order=="o__"~str_c(ncbi.order,"*",sep=""),
                              T~gtdb.order))
write_tsv(order,"0.ncbi2gtdb_order_c.txt")

family = read_tsv("0.ncbi2gtdb_family.txt") %>% 
  mutate(gtdb.family=case_when(gtdb.family=="f__"~str_c(ncbi.family,"*",sep=""),
                               T~gtdb.family))
write_tsv(family,"0.ncbi2gtdb_family_c.txt")

genus = read_tsv("0.ncbi2gtdb_genus.txt") %>% 
  mutate(gtdb.genus=case_when(gtdb.genus=="g__"~str_c(ncbi.genus,"*",sep=""),
                              T~gtdb.genus))
write_tsv(genus,"0.ncbi2gtdb_genus_c.txt")

species = read_tsv("0.ncbi2gtdb_species.txt") %>% 
  mutate(gtdb.species=case_when(gtdb.species=="s__"~str_c(ncbi.species,"*",sep=""),
                                T~gtdb.species))
write_tsv(species,"0.ncbi2gtdb_species_c.txt")

# using GTDBvsNCBI, correct manually
binInfo = read_tsv("1.Bin409.stats.tsv") %>%
  left_join(read_tsv("0.gtdb2ncbi_phylum.txt")) %>%
  left_join(read_tsv("0.gtdb2ncbi_class.txt")) %>%
  left_join(read_tsv("0.gtdb2ncbi_order.txt")) %>%
  left_join(read_tsv("0.gtdb2ncbi_family.txt")) %>%
  left_join(read_tsv("0.gtdb2ncbi_genus.txt")) %>%
  left_join(read_tsv("0.gtdb2ncbi_species.txt")) %>%
  mutate(ncbi.phylum=case_when(is.na(ncbi.phylum)~"p__",
                               T~ncbi.phylum)) %>%
  mutate(ncbi.class=case_when(is.na(ncbi.class)~"c__",
                              T~ncbi.class)) %>%
  mutate(ncbi.order=case_when(is.na(ncbi.order)~"o__",
                              T~ncbi.order)) %>%
  mutate(ncbi.family=case_when(is.na(ncbi.family)~"f__",
                               T~ncbi.family)) %>%
  mutate(ncbi.genus=case_when(is.na(ncbi.genus)~"g__",
                              T~ncbi.genus)) %>%
  mutate(ncbi.species=case_when(is.na(ncbi.species)~"s__",
                                T~ncbi.species)) %>%
  mutate(ncbi.tax=str_c(gtdb.d,ncbi.phylum,ncbi.class,
                        ncbi.order,ncbi.family,ncbi.genus,
                        ncbi.species,sep=";")) %>%
  select(-c(gtdb.d,gtdb.phylum,gtdb.class,gtdb.order,
            gtdb.family,gtdb.genus,gtdb.species,
            ncbi.phylum,ncbi.class,ncbi.order,
            ncbi.family,ncbi.genus,ncbi.species)) %>%
  unique()
# One 2 one
gtdbmap = read_tsv("bac120_metadata_r95.tsv")

# using NCBIvsGTDB, not good
binInfo = read_tsv("1.Bin409.stats.tsv") %>%
  left_join(read_tsv("0.ncbi2gtdb_phylum.txt")) %>%
  left_join(read_tsv("0.ncbi2gtdb_class.txt")) %>%
  left_join(read_tsv("0.ncbi2gtdb_order.txt")) %>%
  left_join(read_tsv("0.ncbi2gtdb_family.txt")) %>%
  left_join(read_tsv("0.ncbi2gtdb_genus.txt")) %>%
  left_join(read_tsv("0.ncbi2gtdb_species.txt")) %>%
  mutate(ncbi.phylum=case_when(is.na(ncbi.phylum)~"p__",
                               T~ncbi.phylum)) %>%
  mutate(ncbi.class=case_when(is.na(ncbi.class)~"c__",
                              T~ncbi.class)) %>%
  mutate(ncbi.order=case_when(is.na(ncbi.order)~"o__",
                              T~ncbi.order)) %>%
  mutate(ncbi.family=case_when(is.na(ncbi.family)~"f__",
                               T~ncbi.family)) %>%
  mutate(ncbi.genus=case_when(is.na(ncbi.genus)~"g__",
                              T~ncbi.genus)) %>%
  mutate(ncbi.species=case_when(is.na(ncbi.species)~"s__",
                                T~ncbi.species)) %>%
  mutate(ncbi.tax=str_c(gtdb.d,ncbi.phylum,ncbi.class,
                        ncbi.order,ncbi.family,ncbi.genus,
                        ncbi.species,sep=";")) %>%
  select(-c(gtdb.d,gtdb.phylum,gtdb.class,gtdb.order,
            gtdb.family,gtdb.genus,gtdb.species,
            ncbi.phylum,ncbi.class,ncbi.order,
            ncbi.family,ncbi.genus,ncbi.species)) %>%
  unique()
write_tsv(binInfo,"2.Bin409.stats.tsv")
