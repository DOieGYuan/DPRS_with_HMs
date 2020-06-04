megahit --k-list 25,31,41,51,61,71,81,101,121,141 -1 *_1.fq.gz -2 *_2.fq.gz -o assembly --out-prefix Coasm -m 0.95 -t 64
# estimate assembly quality
quast -o quast --min-contig 200 --threads 64 Coasm.fa
