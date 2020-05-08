for f in *_1.fq.gz; do kraken2 --db NCBI_Sludge_db/ --paired $f ${f%_1.fq.gz}_2.fq.gz --threads 64 --report ${f%_1.fq.gz}.txt --use-mpa-style; done





