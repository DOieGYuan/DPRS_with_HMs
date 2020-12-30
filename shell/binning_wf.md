# Binning work follow
Here, we explain the binning workflow step-by-Step  
First, we have to activate the conda environment *binning*
```
conda activate binning
# Redirct to work directory
test -e Bining || mkdir Binning && cd Binning
# Link data to the working directory
ln -s ../assembly/Coasm.contigs.fa .
ln -s [reads directory]/*.fq.gz .
```
## Build .SAM and .BAMs files for binning_wf
Use bowtie2 to map metagenomic reads back to the Assembly  
1. Calculate the number of contigs in Coasm.contigs.fa
```
nuc_num=`egrep -v ">" Coasm.contigs.fa | wc -c`
```
2. If $nuc_num > 4000000000, then --large-index is needed
```
if (($nuc_num >= 4000000000)); then
	bowtie2-build Coasm.contigs.fa bz2/asm --threads 64 --large-index
else
	bowtie2-build Coasm.contigs.fa bz2/asm --threads 64
fi
```
3. Align using bowtie2
```
for f in *_1.fq.gz
do bowtie2 -1 $f -2 ${f%_1.fq.gz}_2.fq.gz -x bz2/asm -S ${f%_1.fq.gz}.sam -p 64
done
```
4. Convert .SAM to sorted .BAM using samtools
```
samtools faidx Coasm.contigs.fa
for i in *.sam
do
   samtools import Coasm.contigs.fa $i ${i%.sam}.bam
   samtools sort ${i%.sam}.bam -o ${i%.sam}.sorted.bam -@ 64 -m 4G
   samtools index ${i%.sam}.sorted.bam -@ 64
done
```
## Binning using metabat2
1. Calculate the depth file
```
jgi_summarize_bam_contig_depths --outputDepth depth_metabat.txt *.sorted.bam
```
2. Binning
```
metabat2 -i Coasm.contigs.fa -a depth_metabat.txt -o metabat_binning/metabat -m 1500 -t 64
```
## Binning using maxbin2
1. Use the depth file, but need re-format it to the maxbin2-recognizable style
```
sample_num=`ls *_1.fq.gz | wc -l`
i=1
sum=0
while ((i <= sample_num))
	do ((sum = 2*i+2))
	cut -f 1,$sum depth_metabat.txt | sed '1d' > A${i%}.coverage.tab
	((i++))
done
ls *coverage.tab > abundance.list
```
2. Binning
```
mkdir maxbin
run_MaxBin.pl -contig Coasm.contigs.fa -abund_list abundance.list -out maxbin/maxbin -thread 64 -min_contig_length 1000 -max_iteration 55
mkdir maxbin2_bins
mv maxbin/maxbin.*.fasta maxbin2_bins
```
## Binning using concoct
1. Cut up the fasta
```
cut_up_fasta.py Coasm.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
```
2. Calculate coverage table and PCA
```
concoct_coverage_table.py contigs_10K.bed *.sorted.bam > coverage_table_concoct.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table_concoct.tsv -b concoct -t 64
```
3. Binning
```
mkdir concoct_output
merge_cutup_clustering.py concoct_clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/concoct_bins
extract_fasta_bins.py Coasm.contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/concoct_bins/
```
Clean the workspace
```
rm *.sam
#rm *.bam*
rm concoct_*
rm contigs_10K*
rm abundance.list
rm A*.coverage.tab
#rm *.fai
rm -rf bz2
rm -rf maxbin2
rm -rf concoct
```
## Aggregate the three sets of bins
```
metawrap bin_refinement -o metawrap -t 64 -m 250 -c 50 -x 10 -A maxbin2_bins -B metabat_binning -C concoct_bins
mv metawrap/metawrap_*bins .
rm -rf metawrap
mv metawrap_*bins metawrap/
```
## Refine the bins using refineM
Remove outliers in each MAG (contigs had GC or tetranucleotide distance outside the 98th percentile of the expected distributions or conflicting phylum-level taxonomy with specific MAGs)
1. Removing contamination based on genomic properties
```
test -e refineM_workfile || mkdir refineM_workfile
refinem scaffold_stats -c 64 -x fa Coasm.contigs.fa metawrap refineM_workfile/scaffold_stats *.sorted.bam
refinem outliers refineM_workfile/scaffold_stats/scaffold_stats.tsv refineM_workfile/outlier --cov_corr 0.8 --cov_perc $((nuc_num*50)) # 0.5*sampleNUM*100
refinem filter_bins metawrap refineM_workfile/outlier/outliers.tsv refineM_workfile/filtered_output_params -x fa
```
2. Removing contamination based on taxonomic assignments
```
refinem call_genes -c 64 refineM_workfile/filtered_output_params refineM_workfile/gene_output -x fa
refinem taxon_profile -c 64 refineM_workfile/gene_output refineM_workfile/scaffold_stats/scaffold_stats.tsv $REFINEM_DATABASE/genome_db.2020-07-30.genes.faa [$REFINEM_DATABASE]/gtdb_r95_taxonomy.2020-07-30.tsv  refineM_workfile/taxon_profile
refinem taxon_filter -c 64 refineM_workfile/taxon_profile refineM_workfile/taxon_filter.tsv
refinem filter_bins filtered_output refineM_workfile/taxon_filter.tsv refineM_bins -x fa

# Please manually filter potential incongruent 16S rRNA genes
# Removing contigs with incongruent 16S rRNA genes
refinem ssu_erroneous -c $2 -x fa refineM_bins refineM_workfile/taxon_profile $REFINEM_DATABASE/gtdb_r80_ssu_db.2018-01-18.fna $REFINEM_DATABASE/gtdb_r80_taxonomy.2017-12-15.tsv ssu_output
rm -rf refineM_workfile
```
## reassemble using metaWRAP
1. Combine all the samples
```
zcat *_1.fq.gz > comb_1.fastq
zcat *_2.fq.gz > comb_2.fastq
```
2. Reassemble
```
# if need, change numbers behind -c and -x for controlling the MAGs' quality
metawrap reassemble_bins -b refineM_bins -o reasm -1 comb_1.fastq -2 comb_2.fastq -c 50 -x 10 -t $2 -m 250 -l 500
conda deactivate
rm comb_?.fastq
rm Coasm.contigs.fa
rm -rf reasm/work_files
rm -rf reasm/reassembled_bins.checkm
rm -rf reasm/original_bins
```
Now, binning workflow all done!
