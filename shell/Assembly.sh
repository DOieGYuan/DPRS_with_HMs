forwardseq=$(ls *_1.fq.gz | sed ':a;N;s/\n/,/g;ta')
reverseseq=$(ls *_2.fq.gz | sed ':a;N;s/\n/,/g;ta')
megahit --k-list 27,33,41,51,61,71,81,91,101,111,121,141 --min-contig-len 1000 -1 $forwardseq -2 $reverseseq -o assembly --out-prefix Coasm -m 0.95 -t 64
# Simplified headers
cat assembly/Coasm.contigs.fa | sed 's/ .*//g' > Coasm.contigs.fa
# estimate assembly quality
quast -o quast --min-contig 200 --threads 64 Coasm.fa
