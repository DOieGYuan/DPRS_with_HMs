./download_eggnog_data.py --data_dir db -y -f bact arch
./emapper.py -i data/proteins.faa -o data -d bact --data_dir db/ --usemem --cpu 64 -m diamond

