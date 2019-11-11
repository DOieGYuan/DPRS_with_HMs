test -e hmmsearch_result/ || mkdir hmmsearch_result
for gene in *.hmm;
do
for f in *.faa; do ls $f > hmmsearch_result/${f%.faa}_tmp1; hmmsearch -E 0.00001 --domE 0.00001 --cpu 12 $gene $f | grep domZ >> hmmsearch_result/${f%.faa}_tmp1; cat hmmsearch_result/${f%.faa}_tmp1 | sed ':a;N;s/\.faa\n//g;ta' | sed 's/Domain search space  (domZ)\: */\t/g' | sed 's/ *\[.*//g' > hmmsearch_result/${f%.faa}_tmp; rm hmmsearch_result/${f%.faa}_tmp1; done
cat hmmsearch_result/*_tmp > hmmsearch_result/${gene%.hmm}.tsv
rm hmmsearch_result/*_tmp
done
./merge.py hmmsearch_result/*.tsv > hmmsearch_result/result.txt
sed -i 's/bin\./MAG-/g' hmmsearch_result/result.txt