# Reproduce the results in our paper [Unpublished]
This study focuses on the taxon-specifc heavy metal resistance in denitrifying phosphorus removal sludge
## Main contents in this repository
Home-made functional genes database is deposited in [Database](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Database) filefold (both .fasta, .dmnd and .hmm);  
Raw data are in folder [Raw_data](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Rawdata);  
R scripts to reproduce the figures in our paper are in [Rscripts](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/Rscripts);  
Shell scripts to reproduce data processing in our paper are in [Shellscripts](https://github.com/DOieGYuan/DPRS_with_HMs/tree/master/shell).  

## Install dependencies
All the processes were performed on ubuntu 16.04LTS OS.
We recommond using [Anaconda](https://www.anaconda.com/) to install all the dependencies.   
At least 256GB RAM for assembly, 132GB RAM for downstream taxonomic analysis, and 64GB RAM for the rest analysis.
### quality control
```
conda install -c bioconda fastqc trimmomatic
```
### Assembly
```
conda install -c bioconda megahit quast
```
### Binning
```
conda install -c bioconda metabat2 maxbin2 concoct checkm-genome bowtie2 drep
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
conda install -c bioconda gtdbtk singlem kraken2
```
### Phylogenetic tree construction
Use [UBCG](https://www.ezbiocloud.net/) for building tree file and [iTOL](https://itol.embl.de/) for visualization.
### Functional annotation
```
conda install -c bioconda enrichm hmmer diamond
```
Also install [emapper](https://github.com/eggnogdb/eggnog-mapper) and [InterProScan](https://github.com/ebi-pf-team/interproscan) manually
### R packages
* DESeq2
* vegan
* tidyverse
* ggpubr
* igraph
### Enrichment analysis
We use [GSEA](https://www.gsea-msigdb.org/gsea/) for enrichment analysis and identify the leading-edge subset.  
[Cytoscape](https://cytoscape.org/) and [EnrichmentMap App](http://apps.cytoscape.org/apps/enrichmentmap) are used for visualization following this awesome [protocol](https://doi.org/10.1038/s41596-018-0103-9).
### Co-occurence network
We recommond using [CoNet](http://psbweb05.psb.ugent.be/conet/microbialnetworks/conet_new.php) for predict both the co-presence and mutal-exclusion patterns.
### Metaproteomic analysis
Since the commercial software DA2 maybe inaccessible to some readers, we recommond some other outstanding open-source software for MS data processing:
* [Maxquant](http://www.coxdocs.org/doku.php?id=maxquant:start)
* [COSS](https://github.com/compomics/COSS)
* [MPA](https://github.com/compomics/meta-proteome-analyzer)
* [Peptide-shaker](https://github.com/compomics/peptide-shaker)(for interpretation)  

We encourage readers to use aboved software to process our data to see whether different approaches can produce the same results. Also, the results generated by DA2 are given in [here](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/Rawdata/Metaproteome/Annotation_enrichM.tsv).  
## Work flow
### Download raw data
see **Data availability** section in our paper to download raw data.  
Note that this work flow can also be applied to other data by renaming the paired-end sequences in the form of `[sample_name1]_1.fq.gz` and `[sample_name1]_2.fq.gz`
### Quality control
```
cd [reads folder in your work station]
for f in *.fq.gz
do fastqc $f
done
pip install multiqc
multiqc .
```
To understand the results generated by fastqc, see [Fastqc website](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [video tutorial](http://www.youtube.com/watch?v=bz93ReOv87Y).  
Then, remove adapters and low-quality reads
```
curl -O -L http://dib-training.ucdavis.edu.s3.amazonaws.com/mRNAseq-semi-2015-03-04/TruSeq2-PE.fa
for f in *_1.fq.gz
do TrimmomaticPE $f \
    ${f%_1.fq.gz}_2.fq.gz \
    ${f%_1.fq.gz}_1.qc.fq.gz s1_se \
    ${f%_1.fq.gz}_2.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25
done
```
You can skip this QC process by downloading the clean version of our data (see **Data availability** section in our paper).
### Assembly
Copy [Assembly.sh](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/shell/Assembly.sh) to reads directory and run `./Assembly.sh` (you may need run `chmod 775 Assembly.sh` first) or
```
mkdir assembly
cd assembly
ln -s [reads directory]/*.fq.gz .
megahit --k-list 25,31,41,51,61,71,81,101,121,141 -1 *_1.fq.gz -2 *_2.fq.gz -o assembly --out-prefix Coasm -m 0.95 -t 64
# estimate assembly quality
quast -o quast --min-contig 200 --threads 64 Coasm.contigs.fa
```
Now we obtain **[Coasm.contigs.fa](https://submit.ncbi.nlm.nih.gov/subs/wgs_batch/SUB6248023)**.
### Binning
Use our packaged [binning_wf.sh](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/shell/binning_wf.sh)
```
mkdir binning
cd binning
ln -s ../assembly/Coasm.fa .
ln -s [reads directory]/*.fq.gz .
chmod 775 binning_wf
./binning_wf.sh Coasm.fa [threads] # in our case, 64
```  
Now we get **high-quality metagenomic-assembled genomes (MAGs)**.  
Or, you can skip this time-consuming step by downloading MAGs from [here](https://submit.ncbi.nlm.nih.gov/subs/wgs_batch/SUB6626062/overview).
### Taonomic classification
See [GTDBtk manual](https://ecogenomics.github.io/GTDBTk/) for details in MAG taxonomic classification
Here we use:
```
mkdir genomes
cp -r metawrap/reasm/reassembled_bins/ genomes/
# First, download database
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz
tar xvzf gtdbtk_r89_data.tar.gz
# Then, annotate
gtdbtk classify_wf --genome_dir genomes --out_dir GTDBtk_tax --cpus 64
```
Now we have the taxonomic information for each MAG.
Then we build Phylogenetic tree for MAGs using UBCG pipeline
```
mkdir UBCG
cd ubcg
ln -s metawrap/reasm/reassembled_bins/*.fa .
#Step 1: Converting genome assemblies or contigs (fasta) to bcg files
for f in *.fa ; do java -jar UBCG.jar extract -i $f -bcg_dir ./bcg -label ${f%.fna} -t 64; done
#Step 2: Generating multiple alignments from bcg files
java -jar UBCG.jar align -bcg_dir ./bcg -out_dir tree/ -t 64 -prefix dprs -raxml
cat tree/dprs*92*.nwk | sed 's/\'//g' > TreeFile.nwk
#calculate orthoANI
java -jar OAU.jar -u usearch -fd ./ -n 64 -fmt matrix -o ANI.txt
```
Now we get **[TreeFile.nwk](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/Rawdata/Metagenome/iTOL_PhylogeneticTree/TreeFile.nwk)**  
Upload the .nwk file onto the iTOL online system to construct a phylogenetic tree.  

*(Optional step for Supplementary information Fig.S1)* Short-read based classification using [kraken2](https://ccb.jhu.edu/software/kraken2/)  
Database is based on all the bacterial, algal and viral genomes download from NCBI. (Custom database construction see [kraken2's manual](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual))
```
for f in *_1.fq.gz
do kraken2 --db NCBI_db/ --paired $f ${f%_1.fq.gz}_2.fq.gz \
  --threads 64 --report ${f%_1.fq.gz}.txt --use-mpa-style
done
```
*(Optional step for Supplementary information Fig.S1)* 16S rRNA DNA-based classification using [singleM](https://github.com/wwood/singlem) using [singleM.sh](https://raw.githubusercontent.com/DOieGYuan/DPRS_with_HMs/master/shell/SingleM.sh), we can get 16S-based profile of taxonomy in DPRS and the proportion of community covered by our MAGs.


### Extract functional genes in MAGs
Download our home-made referential database (.dmnd, hmm and original .fasta).  
```
mv [directory]/*.dmnd .
mv [directory]/*.hmm .
# predict CDS by Prodigal
for f in genomes/*.fa
do prodigal -p single -i $f -a ${f%.fa}.faa
done
./Search_functional_genes_DIAMOND.sh # 1e-10 50 80 80
./Search_functional_genes_HMMER.sh 60
```
Couple the results of DIAMOND (1e-10,id50,cov80) and HMMER (socre > 60), the presence of functional genes in MAGs is profiled.  
According to the [annotation template](https://itol.embl.de/help/dataset_binary_template.txt) provided by iTOL, we get **[dataset_binary_functional_genes.txt](https://github.com/DOieGYuan/DPRS_with_HMs/raw/master/Rawdata/Metagenome/iTOL_PhylogeneticTree/dataset_binary_functional_genes.txt)**  

Draw the file to the iTOL tree decorating window to activate this annotation.

### Estimate the abundance of each MAG based on its coverage
We use the "quant_bins" module of metaWRAP to get copies of genome per million reads (CoPM):  
```

```
Then perform LEfSe analysis to identify biomarkers in each group
```

```
