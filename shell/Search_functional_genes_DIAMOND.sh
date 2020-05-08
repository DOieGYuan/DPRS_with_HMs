test -e Diamond_funtion_results||mkdir Diamond_funtion_results
read -p "e-value: " evalue
read -p "identity: " identity
read -p "qcover: " qcover
read -p "scover: " scover

for db in *.dmnd;
do
for f in *.faa; 
do diamond blastp --threads 64 --db ${db%.dmnd} -f 6 qseqid sseqid evalue bitscore length pident nident -q $f --top 3 -e $evalue --id $identity --query-cover $qcover --subject-cover $scover --sensitive | cut -d "	" -f 7 > /dev/shm/${f%.faa}.txt; test -s /dev/shm/${f%.faa}.txt && ls $f | sed 's/bin\./MAG-/g' | sed 's/\.faa//g' >> Diamond_funtion_results/${db%.dmnd}.txt; rm /dev/shm/${f%.faa}.txt; done
sed -i 's/$/\tYes/g' Diamond_funtion_results/${db%.dmnd}.txt
done
./merge.py Diamond_funtion_results/*.txt > Diamond_funtion_results/diamond_results_combined.txt
