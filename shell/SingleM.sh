for f in *_1.fq.gz; do singlem pipe --forward $f --reverse ${f%_1.fq.gz}_2.fq.gz --otu_table ${f%_1.fq.gz}.txt --threads 64; done

# refer to github user manual for details of downstream analysis

singlem summarise --input_otu_tables *.txt --output_otu_table combined_otu.txt
singlem summarise --input_otu_tables combined_otu.txt --krona krona.html
singlem summarise --input_otu_tables s*.txt --cluster --clustered_output_otu_table clustered_otu.txt
singlem summarise --input_otu_tables s*.txt --rarefied_output_otu_table rarefied_otu.txt --number_to_choose 100

# determine how much of a community is represented in an assembly or represented by a set of genomes
singlem pipe --sequences /media/linyuan/SSD1/Binning/metaWarp_refine/reasm/reassembled_bins/*.fa --otu_table genomes_otu_table.txt
singlem pipe --sequences /media/linyuan/SSD1/Binning/all_coasm_kmin21.contigs.fa --otu_table assembly_otu.txt
singlem appraise --metagenome_otu_tables combined_otu.txt --assembly_otu_tables assembly_otu.txt --genome_otu_tables genomes_otu_table.txt
