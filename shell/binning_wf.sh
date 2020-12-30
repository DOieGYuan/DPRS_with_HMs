# usage: source binning_wf.sh [assembly xxx.fasta or xxx.fa] [threads]
# author DOieGYuan
# Github: http://github/DOieGyuan/DPRS_with_HMs/
echo "now building assembly index"
test -e bz2/ || mkdir bz2
#align using bowtie2
nuc_num=`egrep -v ">" $1 | wc -c`
sample_num=`ls *_1.fq.gz | wc -l`
if (($nuc_num >= 4000000000)); then
	bowtie2-build $1 bz2/asm --threads $2 --large-index
else
	bowtie2-build $1 bz2/asm --threads $2
fi
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
# maxbin2
echo "binning by maxbin2"
i=1
sum=0
while ((i <= sample_num))
	do ((sum = 2*i+2))
	cut -f 1,$sum depth_metabat.txt | sed '1d' > A${i%}.coverage.tab
	((i++))
done
ls *coverage.tab > abundance.list
mkdir maxbin
run_MaxBin.pl -contig $1 -abund_list abundance.list -out maxbin/maxbin -thread $2 -min_contig_length 1000 -max_iteration 55
mkdir maxbin2_bins
mv maxbin/maxbin.*.fasta maxbin2_bins
# clean workplace
rm *.sam
rm concoct_*
rm contigs_10K*
rm abundance.list
rm A*.coverage.tab
rm -rf bz2
rm -rf maxbin
mv concoct_output/concoct_bins .
rm -rf concoct_output
conda activate metawrap
# metaWrap
echo "use metaWRAP to refine bins"
metawrap bin_refinement -o metawrap -t $2 -m 120 -c 50 -x 10 -A maxbin2_bins -B metabat_binning -C concoct_bins
mv metawrap/metawrap_*bins .
rm -rf metawrap
mv metawrap_*bins metawrap/
# refineM
source refineM.sh ${1%}.contigs.fa $2 $1
rm *.bam
rm *.fai
# reassemble using metaWRAP
zcat *_1.fq.gz > comb_1.fastq
zcat *_2.fq.gz > comb_2.fastq
# if need, change numbers behind -c and -x for controlling the MAGs' quality
metawrap reassemble_bins -b refineM_bins -o reasm -1 comb_1.fastq -2 comb_2.fastq -c 50 -x 10 -t $2 -m 250 -l 500
conda deactivate
rm comb_?.fastq
rm ${1%}.contigs.fa
rm -rf reasm/work_files
rm -rf reasm/reassembled_bins.checkm
rm -rf reasm/original_bins
echo "binning workflow finished successfully"
