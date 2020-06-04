# usage: ./sh [assembly xxx.fasta or xxx.fa] [threads]
# author DOieGYuan
# Github: http://github/DOieGyuan/DPRS_with_HMs/
echo "now building assembly index"
test -e bz2/ || mkdir bz2
#align using bowtie2
bowtie2-build $1 bz2/asm --threads $2 
echo "align using bowtie2"
for f in *_1.fq.gz
do bowtie2 -1 $f -2 ${f%_1.fq.gz}_2.fq.gz -x bz2/asm -S ${f%_1.fq.gz}.sam -p $2
done
#SAM to BAM using samtools
echo "parsing .bam files"
samtools faidx $1
for i in *.sam
do
   samtools import $1 $i ${i%.sam}.bam
   samtools sort ${i%.sam}.bam -o ${i%.sam}.sorted.bam -@ $2 -m 4G
   samtools index ${i%.sam}.sorted.bam -@ $2
done
# metabat2
echo "binning by metabat2"
jgi_summarize_bam_contig_depths --outputDepth depth_metabat.txt *.sorted.bam
metabat2 -i $1 -a depth_metabat.txt -o metabat_binning/metabat -m 1500 -t $2
# concoct
echo "binning by concoct"
cut_up_fasta.py $1 -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed *.sorted.bam > coverage_table_concoct.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table_concoct.tsv -b concoct -t $2
mkdir concoct_output
merge_cutup_clustering.py concoct_clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/concoct_bins
extract_fasta_bins.py $1 concoct_output/clustering_merged.csv --output_path concoct_output/concoct_bins/
# maxbin2 (change below as needed)
echo "binning by maxbin2"
cut -f 1,4 depth_metabat.txt | sed '1d' > A1.coverage.tab
cut -f 1,6 depth_metabat.txt | sed '1d' > A2.coverage.tab
cut -f 1,8 depth_metabat.txt | sed '1d' > A3.coverage.tab
cut -f 1,10 depth_metabat.txt | sed '1d' > A4.coverage.tab
#cut -f 1,12 depth_metabat.txt | sed '1d' > A5.coverage.tab
#cut -f 1,14 depth_metabat.txt | sed '1d' > A6.coverage.tab
#cut -f 1,16 depth_metabat.txt | sed '1d' > A7.coverage.tab
#cut -f 1,18 depth_metabat.txt | sed '1d' > A8.coverage.tab
#cut -f 1,20 depth_metabat.txt | sed '1d' > A9.coverage.tab
#cut -f 1,22 depth_metabat.txt | sed '1d' > A10.coverage.tab
#cut -f 1,24 depth_metabat.txt | sed '1d' > A11.coverage.tab
#cut -f 1,26 depth_metabat.txt | sed '1d' > A12.coverage.tab

ls *coverage.tab > abundance.list
mkdir maxbin
run_MaxBin.pl -contig $1 -abund_list abundance.list -out maxbin/maxbin -thread $2 -min_contig_length 1000 -max_iteration 55
