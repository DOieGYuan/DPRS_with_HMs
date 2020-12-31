# Protein identification
Downlaod NR database from NCBI
1. Build .dmnd
```
diamond makedb --db nr --in nr.faa --threads 64
```
2. Search homologue sequences
```
diamond blastp --db nr -q Coasm.faa -f 6 sseqid --id 95 --top 3 --scover 80 --qcover 80 -e 1e-10 > homologue.txt
./getID_from_fasta.py -i homologue -f nr.faa -o homologue.faa
```
3. Combine predicted ORFs and homologue sequences
```
cat Coasm.faa homologue.faa > Metap_reference.faa
```
4. Search in PD2.2  
5. We get [Protein match results](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Data/Metaproteome/PD2.2_Match_Results.tsv) and [Peptide match results](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Data/Metaproteome/PeptideGroups.tsv). Besides, the .faa file for identified entries is [metaproteome.faa](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Data/Metaproteome/Metaproteome.faa).  
6. Extract Peptide information 
```
library(tidyverse)

# Extract peptides and proteins linkages
pepinfo = read_tsv("PeptideGroups.txt")
peptide2protein = pepinfo %>%
  dplyr::select(Peptide = peptide, Protein = protein) %>%
  na.omit() %>%
  separate(Protein, paste0("X", 1:(max(str_count(.$Protein,"; "))+1), seq = ""), sep = "; ") %>%
  gather(key = "X", value = "Protein", -Peptide) %>%
  dplyr::select(Peptide, Protein) %>%
  na.omit() %>%
  base::unique()
write_tsv(peptide2protein %>%
            arrange(Protein) %>%
            filter(Protein!=""),
          "Peptide2protein.txt")
```
