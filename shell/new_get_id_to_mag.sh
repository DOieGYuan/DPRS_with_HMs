#source script [output_dir]
test -e ${1%}||mkdir ${1%}
read -p "e-value: " evalue
read -p "identity: " identity
read -p "qcover: " qcover
read -p "scover: " scover

for db in *.dmnd;
do
for f in *.faa; 
do diamond blastp --threads 12 --db ${db%.dmnd} -f 6 sseqid -q $f -k 1 -e $evalue --id $identity --query-cover $qcover --subject-cover $scover --sensitive > /dev/shm/${f%.faa}.txt; test -s /dev/shm/${f%.faa}.txt && cat /dev/shm/${f%.faa}.txt | sed ':a;N;s/\n/,/g;ta' > /dev/shm/${f%.faa}.tmp; test -s /dev/shm/${f%.faa}.txt && ls $f | sed 's/bin\./MAG-/g' | sed 's/\.faa//g' >> /dev/shm/${f%.faa}.tmp; rm /dev/shm/${f%.faa}.txt; test -e /dev/shm/${f%.faa}.tmp && sed -i ':a;N;s/\n/\t/g;ta' /dev/shm/${f%.faa}.tmp;
test -e /dev/shm/${f%.faa}.tmp && cat /dev/shm/${f%.faa}.tmp >> ${1%}/${db%.dmnd}.txt; test -e /dev/shm/${f%.faa}.tmp && rm /dev/shm/${f%.faa}.tmp
done
done
cd ${1%}/
for f in *.txt; do cat $f | sed "s/$/\t${f%.txt}/g" >> Gene2MAG.tsv;done
for f in *.txt; do cat $f | cut -f 1 | sed ':a;N;s/\n/,/g;ta' > ${f%.txt}.tmp; cat $f | cut -f 2 | sed ':a;N;s/\n/,/g;ta'>> ${f%.txt}.tmp ; sed -i ':a;N;s/\n/\t/g;ta' ${f%.txt}.tmp; ls $f >> ${f%.txt}.tmp; sed -i ':a;N;s/\n/\t/g;ta' ${f%.txt}.tmp;done
cat *.txt | cut -f 2 | sort -u > uniqueMAG.list
cat *.tmp > Geneset2MAG.tsv
rm *.tmp
cd ../
