#!/bin/bash
# Sophia Poertner
# Run haybaler taxonomy to attribute NCBI taxonomic lineage to Haybaler output files *haybaler.csv in the current dir. See https://github.com/MHH-RCUG/haybaler/
# Usage: bash run_haybaler_tax.sh
# Adjustments for nf_wochenende by Lisa Hollstein


echo "Starting Haybaler taxonomy"


# check if requirements for pytaxonkit are installed 
if ! [ -d /home/$USER/.taxonkit ]; then
  echo "INFO: Haybaler taxonomy: Requirements for pytaxonkit not found. Trying to download/install them now"
  wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -zxvf taxdump.tar.gz && mkdir -p $HOME/.taxonkit && cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit || echo "failed to install requirements" && exit
fi


# for samples
count=`ls -1 *haybaler.csv 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then
    for csv in $(ls *haybaler.csv)
    do
      python3 haybaler_taxonomy.py -i $csv -p .
    done
fi


# move log files into subdirectory
mkdir -p logs
for file in $(ls no_*_found_for_chromosomes.csv)
do
  mv $file logs/
done


# uncomment the next section for testing references (just for developers!)
count=`ls -1 *fa*.fai 2>/dev/null | wc -l`
if [[ $count != 0 ]]
    then
    for fai in $(ls *fa*.fai)
    do
      python3 haybaler_taxonomy.py -i $fai -p . -t True
    done
fi
##### taxonomy ######
# python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /mnt/ngsnfs/seqres/metagenref/

