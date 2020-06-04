#usage: ./makedmnd.sh [listfile] [db_basename]
cat $1 | sort -u > ${1%}.unique
./getID_from_fasta.py -f metaproteome.fasta -i ${1%}.unique -o ${2%}.faa
diamond makedb --in ${2%}.faa --db ${2%}
rm $1
rm ${1%}.unique
rm ${2%}.faa
mv ${2%}.dmnd found_gene_in_MAG
