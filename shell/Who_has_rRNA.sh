printf "bin\tSSU_rRNA\n" > SSU_rRNA.txt
for f in *.cmscan.txt;
do cat $f | grep "SSU_rRNA" > /dev/shm/${f%.cmscan.txt}.tmp
test -s /dev/shm/${f%.cmscan.txt}.tmp && ls $f | sed 's/\.cmscan\.txt//g' | sed 's/$/\tPresent/g' >> SSU_rRNA.txt
done
rm /dev/shm/*.tmp

printf "bin\tLSU_rRNA\n" > LSU_rRNA.txt
for f in *.cmscan.txt;
do cat $f | grep "LSU_rRNA" > /dev/shm/${f%.cmscan.txt}.tmp
test -s /dev/shm/${f%.cmscan.txt}.tmp && ls $f | sed 's/\.cmscan\.txt//g' | sed 's/$/\tPresent/g' >> LSU_rRNA.txt
done
rm /dev/shm/*.tmp

printf "bin\t5S_rRNA\n" > 5S_rRNA.txt
for f in *.cmscan.txt;
do cat $f | grep "5S_rRNA" > /dev/shm/${f%.cmscan.txt}.tmp
test -s /dev/shm/${f%.cmscan.txt}.tmp && ls $f | sed 's/\.cmscan\.txt//g' | sed 's/$/\tPresent/g' >> 5S_rRNA.txt
done
rm /dev/shm/*.tmp
