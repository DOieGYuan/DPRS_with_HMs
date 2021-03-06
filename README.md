# Revealing taxon-specific heavy metal resistance mechanisms in denitrifying phosphorus removal sludge using genome-centric metaproteomics [*Microbiome*](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01016-x)
Follow this step-by-step guidance (from raw data processing to plotting) to reproduce the bioinformatic results in our paper.  
The background, significance, and main contents (findings) in this study is extensively disscussed in the paper.  

## Main contents in this repository
Functional genes references (.fasta and .dmnd files) are deposited in [Database](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Database);  
Raw data, intermediate files (if applicable) and final results are in folder [Data](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Data);  
R scripts to reproduce the figures in our paper are in [Rscripts](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Rscripts);  
Shell scripts to reproduce data processing in our paper are in [Shellscripts](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/shell).  
* **[Install dependencies](https://github.com/DOieGYuan/DPRS_with_HMs#install-dependencies)**
* [quality control](https://github.com/DOieGYuan/DPRS_with_HMs#quality-control)
* [Assembly](https://github.com/DOieGYuan/DPRS_with_HMs#assembly)
* [Binning](https://github.com/DOieGYuan/DPRS_with_HMs#binning)
* [Taxonomic classification](https://github.com/DOieGYuan/DPRS_with_HMs#taxonomic-classification)
* [Phylogenetic tree construction](https://github.com/DOieGYuan/DPRS_with_HMs#phylogenetic-tree-construction)
* [Functional annotation](https://github.com/DOieGYuan/DPRS_with_HMs#functional-annotation)
* [R packages](https://github.com/DOieGYuan/DPRS_with_HMs#r-packages)
* [Enrichment analysis](https://github.com/DOieGYuan/DPRS_with_HMs#enrichment-analysis)
* [Co-occurence network](https://github.com/DOieGYuan/DPRS_with_HMs#co-occurence-network)
* [Metaproteomic analysis](https://github.com/DOieGYuan/DPRS_with_HMs#metaproteomic-analysis)
* **[Work flow](https://github.com/DOieGYuan/DPRS_with_HMs#work-flow)**
* [Download raw data](https://github.com/DOieGYuan/DPRS_with_HMs#download-raw-data)
* [Quality control](https://github.com/DOieGYuan/DPRS_with_HMs#quality-control-1)
* [Assembly](https://github.com/DOieGYuan/DPRS_with_HMs#assembly-1)
* [Binning](https://github.com/DOieGYuan/DPRS_with_HMs#binning-1)
* [Taonomic classification](https://github.com/DOieGYuan/DPRS_with_HMs#taonomic-classification)
* [Extract functional genes in MAGs](https://github.com/DOieGYuan/DPRS_with_HMs#extract-functional-genes-in-mags)
* [Estimate the abundance of each MAG based on its coverage](https://github.com/DOieGYuan/DPRS_with_HMs#estimate-the-abundance-of-each-mag-based-on-its-coverage)
* [Plot abundance profile of MAGs](https://github.com/DOieGYuan/DPRS_with_HMs#plot-abundance-profile-of-mags)
* [Functional annoation of MAGs and Metaproteome](https://github.com/DOieGYuan/DPRS_with_HMs#functional-annoation-of-mags-and-metaproteome)
* [Enrichment analysis](https://github.com/DOieGYuan/DPRS_with_HMs#enrichment-analysis-1)
* [Co-occurence network](https://github.com/DOieGYuan/DPRS_with_HMs#co-occurence-network-1)
* **[Contact Us](https://github.com/DOieGYuan/DPRS_with_HMs#contact-us)**

## Install dependencies
All the processes were performed on ubuntu 16.04LTS OS.
We recommond using [Anaconda](https://www.anaconda.com/) to install all the dependencies.   
At least 256GB RAM for assembly, 150GB RAM for downstream taxonomic analysis, and 64GB RAM for the rest analyses.
### quality control
```
conda create -n qc -c bioconda fastqc multiqc trimmomatic -y
```
### Assembly
```
conda create -n assembly install -c bioconda megahit quast -y
```
### Binning
```
conda create -n binning -c bioconda metabat2 maxbin2 concoct checkm-genome bowtie2 refinem diamond krona prodigal blast infernal trnascan-se -y
```
install [metaWRAP](https://github.com/bxlab/metaWRAP)

>conda create -n metawrap python=2.7  
>conda activate metawrap  
>conda config --add channels defaults  
>conda config --add channels conda-forge  
>conda config --add channels bioconda  
>conda config --add channels ursky  
>conda install --only-deps -c ursky metawrap-mg  

### Taxonomic classification
```
conda create tax install -c bioconda gtdbtk singlem
```
Download the latest database (r95) for GTDB-tk
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
tar xvzf gtdbtk_r95_data.tar.gz
```
### Phylogenetic tree construction
Use [UBCG](https://www.ezbiocloud.net/) for building tree file and [iTOL](https://itol.embl.de/) for visualization.
### Functional annotation
```
conda install -c bioconda enrichm hmmer diamond
```
Also install [emapper](https://github.com/eggnogdb/eggnog-mapper), [InterProScan](https://github.com/ebi-pf-team/interproscan), and [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) manually.
### R packages
* DESeq2
* vegan
* tidyverse
* ggpubr
* igraph
### Enrichment analysis
We use [GSEA](https://www.gsea-msigdb.org/gsea/) for enrichment analysis and identify the leading-edge subset (core-functioning microbes).  
[Cytoscape](https://cytoscape.org/) and [EnrichmentMap App](http://apps.cytoscape.org/apps/enrichmentmap) are used for visualization following this awesome [protocol](https://doi.org/10.1038/s41596-018-0103-9).
### Co-occurence network
We recommond using [CoNet](http://psbweb05.psb.ugent.be/conet/microbialnetworks/conet_new.php) for predict both the co-presence and mutal-exclusion patterns.
### Metaproteomic analysis
Since the commercial software PD2 (Thermo Scientific) maybe inaccessible to some readers, we recommond some other outstanding open-source software for MS data processing:
* [Maxquant](http://www.coxdocs.org/doku.php?id=maxquant:start)
* [COSS](https://github.com/compomics/COSS)
* [MPA](https://github.com/compomics/meta-proteome-analyzer)
* [Peptide-shaker](https://github.com/compomics/peptide-shaker)(for interpretation)  
* [Schiebenhoefer et al.](https://doi.org/10.1038/s41596-020-0368-7) provides a detailed and easy-to-follow protocol for classical proteomics workflow.  

We encourage readers to use aboved software to process our data to see whether different approaches can produce the same results. Also, the results generated by PD2 are given [here](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/Rawdata/Metaproteome/PD2.2_Match_Results.tsv).  
## Work flow
### Download raw data
see **Data availability** section in our paper to download raw data (NGS .fq data and MS/MS .raw data).  
Note that this work flow can also be applied to other data by renaming the paired-end sequences in the form of `[sample_name1]_1.fq.gz` and `[sample_name1]_2.fq.gz`  
*[] means you may need to change the content manually accordingly*
### Quality control
```
conda activate qc
cd [reads folder in your work station]
for f in *.fq.gz
do fastqc $f
done
# Combine all the reports
multiqc .
```
To understand the results generated by fastqc, see [Fastqc website](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [video tutorial](http://www.youtube.com/watch?v=bz93ReOv87Y).  
Then, remove adapters and low-quality reads
```
test -e TruSeq3-PE.fa || ln -s ~/anaconda3/envs/qc/share/trimmomatic/adapters/TruSeq3-PE.fa .
for f in *_1.fq.gz;
do echo -e "\033[32mNow trimming ${f%_1.fq.gz}\033[0m"
trimmomatic PE $f ${f%_1.fq.gz}_2.fq.gz -baseout ${f%_1.fq.gz}.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:40:15:1:true SLIDINGWINDOW:4:15 MINLEN:50 -threads 64
mv $f QC/rawReads/${f%_1.fq.gz}.Raw_1.fq.gz
mv ${f%_1.fq.gz}_2.fq.gz QC/rawReads/${f%_1.fq.gz}.Raw_2.fq.gz
mv ${f%_1.fq.gz}_1U.fq.gz QC/unpairedReads/${f%_1.fq.gz}.unpaired_1.fq.gz
mv ${f%_1.fq.gz}_2U.fq.gz QC/unpairedReads/${f%_1.fq.gz}.unpaired_2.fq.gz
mv ${f%_1.fq.gz}_1P.fq.gz ${f%_1.fq.gz}_1.fq.gz
mv ${f%_1.fq.gz}_2P.fq.gz ${f%_1.fq.gz}_2.fq.gz
done
conda deactivate
cd -
```
### Assembly
Copy [Assembly.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/Assembly.sh) to reads directory and run `./Assembly.sh` (you may need run `chmod 775 Assembly.sh` first) or
```
conda activate assembly
mkdir assembly
cd assembly
ln -s [reads directory]/*.fq.gz .
forwardseq=$(ls *_1.fq.gz | sed ':a;N;s/\n/,/g;ta')
reverseseq=$(ls *_2.fq.gz | sed ':a;N;s/\n/,/g;ta')
megahit --k-list 27,33,41,51,61,71,81,91,101,111,121,141 --min-contig-len 1000 -1 $forwardseq -2 $reverseseq -o assembly --out-prefix Coasm -m 0.95 -t 64
# Simplified headers
cat assembly/Coasm.contigs.fa | sed 's/ .*//g' > Coasm.contigs.fa
# estimate assembly quality
quast -o quast --min-contig 1000 --threads 64 Coasm.contigs.fa
conda deactivate
cd -
```
Now we obtain **[Coasm.contigs.fa](https://ndownloader.figshare.com/files/25915446)**.
### Binning
Use our packaged [binning_wf.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/binning_wf.sh)
```
conda activate binning
mkdir binning
cd binning
ln -s ../assembly/Coasm.contigs.fa .
ln -s [reads directory]/*.fq.gz .
chmod 775 binning_wf
./binning_wf.sh Coasm.fa [threads] # in our case, 64
conda deactivate
```
Details for binning see [here](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/binning_wf.md)  
Now we get **MIMAG high- or medium-quality metagenomic-assembled genomes (MAGs)**.  
Or, skip this time-consuming step by downloading MAGs (entitled "new MAG-xxx") from [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA592128).
### Taonomic classification
See [GTDBtk manual](https://ecogenomics.github.io/GTDBTk/) for details in MAG taxonomic classification  
Here we use:
```
conda activate gtdb
mkdir genomes
cp -r metawrap/reasm/reassembled_bins/ genomes/
gtdbtk classify_wf --genome_dir genomes --out_dir GTDBtk_tax --cpus 64
conda deactivate
```
Now we have the GTDB taxonomic information for each MAG.  
Convert it into NCBI taxonomy [using gtdb_vs_ncbi_r95_bacteria.xlsx](https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdb_vs_ncbi_r95_bacteria.xlsx) using [Convert_GTDB2NCBI.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/Convert_GTDB2NCBI.R). **Manual verification is needed**
### Determine the quality of MAGs
Estimate the completeness and contamination
```
conda activate binning
cd [directory with MAGs]
checkm lineage_wf -t 12 -x fa ./ ./checkm_qa > bin_quality.log
```
Estimate the presence of tRNAs and rRNAs using [How_many_tRNA.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/How_many_tRNA.sh) and [Who_has_rRNA.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/Who_has_rRNA.sh), or manually,
```
# scan rRNA
for f in *.fa;do cmscan --cpu 12 --rfam --cut_ga --nohmmonly --tblout ../cmscan_results/${f%.fa}.cmscan.txt --fmt 2 --clanin ../Rfam.clanin /media/linyuan/HD_Files/Database/Rfam/Rfam.cm $f;done

# scan tRNA
for f in *.fa;do tRNAscan-SE -B --brief --thread 12 -o ../trnascan_results/${f%.fa}.trna.txt $f;done
```
Based on the above information, MAGs are catergoried into high-quality (completeness > 90%, contamination < 5%, presence of the 23S, 16S, 5S rRNA genes, and at least 18 tRNAs), medium-quality (completeness ≥ 50% and contamination < 10%), and low-quality (the rest) based on the [MIMAG standard](https://doi.org/10.1038/nbt.3893). Only MAGs meet with the high- or medium quality and has a quality score ≥ 45 (defined as completeness−5×contamination) were included in the downstream analysis.
### Phylogenetic tree
Build Phylogenetic tree for MAGs using UBCG pipeline
```
conda
mkdir UBCG
cd [ubcg directory]
ln -s [bin_directory]/*.fa .
#Step 1: Converting genome assemblies or contigs (fasta) to bcg files
for f in *.fa ; do java -jar UBCG.jar extract -i $f -bcg_dir ./bcg -label ${f%.fa} -t 64; done
#Step 2: Generating multiple alignments from bcg files
java -jar UBCG.jar align -bcg_dir ./bcg -out_dir tree/ -t 64 -prefix dprs -raxml
cat tree/dprs*92*.nwk | sed 's/\'//g' > TreeFile.nwk
#calculate orthoANI
java -jar OAU.jar -u usearch -fd ./ -n 64 -fmt matrix -o ANI.txt
```
Now we get **[TreeFile.nwk](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Data/iTOL_PhylogeneticTree/TreeFile.nwk)**  
Upload the .nwk file onto the [iTOL online system](https://itol.embl.de/) to construct a phylogenetic tree.  

**Question**: How to decorate the tree to build Fig.1? see [here](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/Decorate_tree.md) for details.  

*(Optional)* Short-read based classification using [kraken2](https://ccb.jhu.edu/software/kraken2/)  
Database is based on all the bacterial, algal and viral genomes download from NCBI. (Custom database construction see [kraken2's manual](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual))
```
for f in *_1.fq.gz
do kraken2 --db NCBI_db/ --paired $f ${f%_1.fq.gz}_2.fq.gz \
  --threads 64 --report ${f%_1.fq.gz}.txt --use-mpa-style
done
```
**(Supplementary information Fig.S1)** ribosomal protein gene-based classification by [singleM](https://github.com/wwood/singlem) using [singleM.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/SingleM.sh), we can get ribosmal protein gene(e.g., *rplB*)-based profile of taxonomy in DPRS and the proportion of community covered by our MAGs.


### Predict functional genes
Predict protein-coding genes (ORFs) using Prodigal,  
1. in co-assembled metagenome (for buildling custom database for protein dentification, see details [here](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/Building_custom_database_for_protein_identification.md))
```
prodigal -a Coasm.faa -d Coasm.fna -i Coasm.fa -p meta
```
2. in MAGs (for exploring the funtional potentials in MAGs)

```
# predict CDS by Prodigal
for f in genomes/*.fa
do prodigal -p single -i $f -a aa/${f%.fa}.faa
done
```
### Annotate genes of interest (those involved in N/P metabolisms and metal resistance) in MAGs
1. DIAMOND searching  
Download our home-made referential database in this repository (.dmnd) and [BacMet](http://bacmet.biomedicine.gu.se/download_temporary.html).
```
cd aa/
mv [directory]/*.dmnd .
./Search_functional_genes_DIAMOND.sh # 1e-10 50 80 80
```
2. Cross validation by METABOLIC based on several Hidden Markov Models (HMMs).
```
conda activate [metabolic]
# Genome files' suffix should be .fasta
perl METABOLIC-G.pl -t 64 -in-gn [genomes] -o metabolic_results -p single
# Extract needed information
cd metabolic_results/Each_HMM_Amino_Acid_Sequence
for f in *.hmm.collection.faa;
do printf "bin\t${f%.hmm.collection.faa}.hmm\n" > ../FunctionalGenesPA/${f%.hmm.collection.faa};
grep ">" $f | sed 's/>//g' | sed 's/~~/\t/g' | cut -f 1 | sort -u | sed 's/$/\tPresent/g' >> ../FunctionalGenesPA/${f%.hmm.collection.faa};
done

# Mergy in R
R
library(tidyverse)
files = dir("FunctionalGenesPA/")
for(f in files){
	nowf = read_tsv(paste("merged/",f,sep=""))
	if(f == files[1]){
		mg = nowf
	}
	else{
		mg=full_join(mg,nowf)
	}
}
write_tsv(mg,"merged/HMM_results.tsv")
```
3. The presence of functional genes in MAGs should be valided by both DIAMOND and METABOLIC.  

### Estimate the abundance of each MAG based on its coverage
We use the "quant_bins" module of metaWRAP to get copies of genome per million reads (CoPM):  
```
unzip *.fq.gz
for f in *.fq
do mv $f ${f%.fq}.fastq
done
metawrap quant_bins -b metawrap/reasm/reassembled_bins/ \
  -o quant_bins/ -a assembly/Coasm.contigs.fa *.fastq -t 64
```
Then base on the [abundance file](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Data/iTOL_PhylogeneticTree/MAG_abundance.tsv), we perform LEfSe analysis to identify biomarkers in each group.  
```
lefse-format_input.py MAG_Abundance.tsv lefse.out -c 2 -s 3 -u 1 -o 1000000
run_lefse.py lefse.out lefse.tsv -l 2
```
### Functional annoation of MAGs and Metaproteome
Please refer to [enrichM](https://github.com/geronimp/enrichM), [emapper](https://github.com/eggnogdb/eggnog-mapper), and [InterProScan](https://github.com/ebi-pf-team/interproscan) to have a more comprehensive understanding of their usages.  
*(Optional) If need KEGG BRITE and KO annotation file, please install R package KEGGREST and use make.BRITE.formatted.R to formate the htext file into tabular .tsv file.*
In our study, we perform the following pipeline:
```
mkdir annotation

# enrichM annotation
enrichm annotate --protein_directory ../genomes/aa --output enrichM \
  --cpus 64 --parallel 64 --ko

# emapper annotation
mkdir emapper/
cd emapper/
download_eggnog_data.py --data_dir db -y -f bact arch
ln -s ../genomes/aa/*.faa .
for f in *.faa;
do python emapper.py -i $f -o ${f%.faa} -d bact -m diamond --cpu 64 --data_dir db/ --usemem
down
cd ..

# interproscan annotation
mkdir interproscan
cd interproscan
ln -s ../genomes/aa/*.faa .
for f in *.faa;
do interproscan.sh -i $f -b ${f%.faa} -cpu 64 -d ./ -f TSV -goterms -hm -iprlookup -pa
done
```
Now, we have metagenome and genomes annotated by enrichM, emapper and interproscan.  
Combine all the annotation into one file and we get the **Additional file 3**.  
**Plot Fig3**: use [plot_Fig.3b&S7.R](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Rscripts/plot_Fig.3b&S7.R) to plot bubble plot to exhibit differentlt expressed functions in nitrifiers (Fig. 3b) and PAOs (Fig. S7)  
### Map protein sequences to MAGs
Use [get_functioal_genes_in_MAGs.sh](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/get_functional_genes_in_MAGs.sh) to get the associations between identified proteins and MAGs.  
```
diamond makedb --in metaprteome.faa --db metaproteome
cd [working directory]
ln -s [bins directory]/*.fa .
./get_genes_in_MAGs.sh Metaproteome
1e-10
95
90
90
```
Then use Peptide2protein.txt build in this [step](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/shell/Building_custom_database_for_protein_identification.md).  
Matches with imperfect peptide alignments (identity < 100 or coverage < 100) should be discarded.  
Now we obtain the information in **Additional file 3**.
### Enrichment analysis
We highly recommond users following the [offical document]((https://www.gsea-msigdb.org/gsea/)) to conduct this step, GSEA has excellent UI and easy to use.  
Here, we provide our study as an example.  
* Use [prepare_GMT.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/Enrichment%20analysis/prepare_GMT.R) to make GMT file and format the expression data as required.
* Pre-rank the expression list.
* Launch the analysis.
* Visualize in Cytoscape, and output the node files and edge files.
* Fill in information (for mapping) for the reloading into cytoscape.
* Make the network clearer.
* Output the network
Now we have **Fig.5**.

### Co-occurence network
Construct the network using CoNet.Please refer to the [manual](http://psbweb05.psb.ugent.be/conet/microbialnetworks/conet_new.php) provided by the author.  
Export the node and edge files.  
Draw Fig. 6ab using [plot_Fig.6ab.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/plot_Fig.6ab.R) and [plot_Fig.6c.R](https://github.com/DOieGYuan/DPRS_with_HMs/blob/master/Rscripts/plot_Fig.6c.R).  
Finally, we finish **Fig.6**.  

Now we have finished the majority of the analyses in our paper.  
*The metagenomic workflow can be finished easily by using our developing [pipeline](https://github.com/DOieGYuan/Easy-genomic-centric-metagenomics-pipepline) with single command*

## Contact Us
Please feel free to contact linyuan@smail.nju.edu.cn  
  
If you find any resource is useful, please cite https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01016-x.  
Lin, Y., Wang, L., Xu, K. et al. Revealing taxon-specific heavy metal-resistance mechanisms in denitrifying phosphorus removal sludge using genome-centric metaproteomics. Microbiome 9, 67 (2021). https://doi.org/10.1186/s40168-021-01016-x
