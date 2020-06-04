# usage: ./script [.dmnd]
test -e ${1%.dmnd}||mkdir ${1%.dmnd}
read -p "e-value: " evalue
read -p "identity: " identity
read -p "qcover: " qcover
read -p "scover: " scover

echo "ID2	contig	nident	genome" > ${1%.dmnd}/${1%.dmnd}.txt
for f in *.faa; 
do diamond blastp --threads 12 --db ${1%.dmnd} -f 6 sseqid qseqid nident full_qseq -q $f -k 1 -e $evalue --id $identity --query-cover $qcover --subject-cover $scover --sensitive > /dev/shm/${f%.faa}.txt; test -s /dev/shm/${f%.faa}.txt && cat /dev/shm/${f%.faa}.txt | sed "s/$/\t${f%.faa}/g" | sed 's/bin\./MAG-/g' > /dev/shm/${f%.faa}.tmp; rm /dev/shm/${f%.faa}.txt;
test -e /dev/shm/${f%.faa}.tmp && cat /dev/shm/${f%.faa}.tmp >> ${1%.dmnd}/${1%.dmnd}.txt; test -e /dev/shm/${f%.faa}.tmp && rm /dev/shm/${f%.faa}.tmp
done
