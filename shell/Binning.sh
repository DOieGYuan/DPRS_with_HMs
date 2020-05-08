#maxbin
#use read_list

run_MaxBin.pl -contig all_coasm_kmin21.contigs.fa -reads_list Coasm_reads -out MaxBin_out -thread 64

#use abundance_list
#abundance files can be achieve from *.sorted.bam

#bedtools

for i in *sorted.bam
do
  genomeCoverageBed -ibam $i > ${i/.pe*/}.histogram.tab
done

#prerequiste
#wget https://raw.githubusercontent.com/ngs-docs/2017-cicese-metagenomics/master/files/calculate-contig-coverage.py
#download pandas

for hist in *histogram.tab
do
  python calculate-contig-coverage.py $hist
done

# or get from depth_metabat.txt
cut -f 1,4 depth_metabat.txt | sed '1d' > A1.coverage.tab
cut -f 1,6 depth_metabat.txt | sed '1d' > A2.coverage.tab
cut -f 1,8 depth_metabat.txt | sed '1d' > A3.coverage.tab
cut -f 1,10 depth_metabat.txt | sed '1d' > A4.coverage.tab
# ...... and so on

ls *coverage.tab > abundance.list
mkdir maxbin
run_MaxBin.pl -contig new_rab_scaffolds.fasta -abund_list abundance.list -out maxbin/maxbin -thread 64 -min_contig_length 1000 -max_iteration 55

#align using bowtie2

bowtie2-build scaffolds.fasta rab --threads 64

for f in *_1.fq.gz
do bowtie2 -1 $f -2 ${f%_1.fq.gz}_2.fq.gz -x rab -S ${f%_1.fq.gz}.sam -p 64
done

#SAM to BAM using samtools

samtools faidx *.fasta

for i in *.sam
do
   samtools import *.fasta $i ${i%.sam}.bam
   samtools sort ${i%.sam}.bam -o ${i%.sam}.sorted.bam -@ 64 -m 8G
   samtools index ${i%.sam}.sorted.bam -@ 64
done


#metabat2

jgi_summarize_bam_contig_depths --outputDepth depth_metabat.txt *.sorted.bam

metabat2 -i scaffolds.fasta -a depth_metabat.txt -o metabat_binning/metabat -m 1500 -t 32

#concoct

cut_up_fasta.py scaffolds.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py contigs_10K.bed *.sorted.bam > coverage_table_concoct.tsv

concoct --composition_file contigs_10K.fa --coverage_file coverage_table_concoct.tsv -b concoct -t 32

mkdir concoct_output

merge_cutup_clustering.py concoct_clustering_gt1000.csv > concoct_output/clustering_merged.csv

mkdir concoct_output/concoct_bins

extract_fasta_bins.py new_rab_scaffolds.fasta concoct_output/clustering_merged.csv --output_path concoct_output/concoct_bins/

#DAS-tool
modify DAS_Tool-1.1.1/src/Fasta_to_Scaffolds2Bin.sh by deleting the first line.


~/App/Genomics/DAS_Tool-1.1.1/src/Fasta_to_Scaffolds2Bin.sh -i maxbin -e fasta > maxbin.scaffolds2bin.tsv
~/App/Genomics/DAS_Tool-1.1.1/src/Fasta_to_Scaffolds2Bin.sh -i metabat -e fa > metabat.scaffolds2bin.tsv
~/App/Genomics/DAS_Tool-1.1.1/src/Fasta_to_Scaffolds2Bin.sh -i concoct -e fa > concoct.scaffolds2bin.tsv

prodigal -p meta -a new_rab_cds.faa -d new_rab_cds.fna -i new_rab_scaffolds.fasta -o new_rab_cds.gbk

DAS_Tool -i maxbin.scaffolds2bin.tsv,concoct.scaffolds2bin.tsv,metabat.scaffolds2bin.tsv -c new_rab_scaffolds.fasta -l maxbin,concoct,metabat -o DASTool/dastool --search_engine diamond --write_bins 1 --proteins new_rab_cds.faa -t 12

#metaWrap
#bin_refine
~/App/Genomics/metaWRAP/bin/metaWRAP bin_refinement -o metawrap_refine -t 64 -A metabat_binning/ -B maxbin2_binning/ -C concoct_binning/ -c 70 -x 10
#vizplot
~/App/Genomics/metaWRAP/bin/metaWRAP blobology -a original_contigs.fa -t 64 -o metaWarp_refine/blobology/ --bins metaWarp_refine/metawrap_70_10_bins *.fastq
#abundance estimate
~/App/Genomics/metaWRAP/bin/metawrap quant_bins -b metaWarp_refine/metawrap_70_10_bins -o metaWarp_refine/quant_bins/ -a original_contigs.fa s*.fastq -t 64
#reassembly
~/App/Genomics/metaWRAP/bin/metaWRAP reassemble_bins -o metawarp/reasm -1 *_1.fastq -2 *_2.fastq -t 64 -m 1280 -c 70 -x 10 -b metawarp/metawrap_70_10_bins
#taxonomy
~/App/Genomics/metaWRAP/bin/metaWRAP classify_bins -b metawarp/reasm/reassembled_bins -o metawarp/bin_tax 
#function_annotation
~/App/Genomics/metaWRAP/bin/metaWRAP annotate_bins -o metawarp/function -t 64 -b metawarp/reasm/reassembled_bins


#checkm

checkm lineage_wf -t 64 ../dastool_DASTool_bins/ ./lineage_wf -x fa
checkm qa ./lineage_wf/lineage.ms ./lineage_wf/ > qa_report.txt
