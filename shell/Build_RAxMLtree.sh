# ./script [fasta] [outdir]
./getID_from_fasta.py -f metaproteome.faa -i ${1%.fasta}.list -o /dev/shm/${1%.fasta}.fasta
muscle -in /dev/shm/$1 -out /dev/shm/${1%.fasta}.ali.fasta
raxmlHPC-PTHREADS -p 12345 -m PROTGAMMALGX -s /dev/shm/${1%.fasta}.ali.fasta -z ${1%.fasta} -n nwk -w [Absolute directory of working file]/${2%} -N 1000 -T 12
rm /dev/shm/${1%.fasta}.ali.fasta
rm /dev/shm/${1%.fasta}.fasta
rm /media/linyuan/SSD1/Annotation/${2%}/RAxML*RUN*
mv /media/linyuan/SSD1/Annotation/${2%}/RAxML_bestTree.nwk [Absolute directory of working file]/${2%}/${1%.fasta}.nwk
mv /media/linyuan/SSD1/Annotation/${2%}/RAxML_info.nwk [Absolute directory of working file]/${2%}/${1%.fasta}_info.txt
