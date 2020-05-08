#Step 1: Converting genome assemblies or contigs (fasta) to bcg files

for f in *.fna ; do java -jar UBCG.jar extract -i $f -bcg_dir ./ -label ${f%.fna} -t 64; done

#Step 2: Generating multiple alignments from bcg files

java -jar UBCG.jar align -bcg_dir ./ -out_dir ubcg/ -t 64 -prefix ubcg -raxml

#calculate orthoANI

conda activate ubcg
java -jar OAU.jar -u usearch -fd dereplicated_genomes/ -n 12 -fmt matrix -o ANI.txt





