conda activate 16S
echo "Which bin?"
read -p "please enter: " genome
read -p "Squences identity: " identity
ID=$(blastn -query sludge_bins/bin.$genome.fa -db SILVA_132_SSUParc_tax_silva.fasta -evalue 0.00001 -outfmt 10 -num_threads 32 -num_alignments 5 -perc_identity $identity -subject_besthit | cut -d "," -f 2 | sed ':a;N;s/\n/\|/g;ta')
cat SILVA_132_SSUParc_tax_silva.fasta | egrep $ID
