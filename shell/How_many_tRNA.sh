printf "bin\ttRNA\n" > tRNA.count.txt
for f in *.tRNA.txt;
do ls $f | sed 's/\.tRNA\.txt//g' > /dev/shm/${f%.tRNA.txt}.tmp
cut -f 5 $f | sed 's/Ile2/Ile/g' | sed 's/fMet/Met/g' | sort -u | egrep -v "Undet" | wc -l >> /dev/shm/${f%.tRNA.txt}.tmp
cat /dev/shm/${f%.tRNA.txt}.tmp | sed ':a;N;s/\n/\t/g;ta' >> tRNA.count.txt
done
rm /dev/shm/bin.*.tmp

