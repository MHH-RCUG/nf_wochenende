#!/bin/bash
echo "Starting the nf_wochenende pipeline"
echo "Requires nextflow in path. In case of problems, activate environment: conda activate wochenende"

#get most current version - for developers
#git pull -q

# ignore system JAVA_HOME, use that supplied by wochenende conda env for nextflow
unset JAVA_HOME
unset WOCHENENDE_DIR
unset HAYBALER_DIR

# cleanup work/ 
rm -rf work

# Developer tests
# get test files, create reference (small files). Do it yourself for bigger test references 
# eg mock community from SRA https://github.com/colindaven/wochenende_manuscript/blob/main/mock/download_fastq.sh
# eg mock ref file human22_zymo_test.fa
#bwa index test/data/ref.fa
#cp -f test/data/*.fastq . && bwa index test/data/ref.fa
#cp -f test/data/*.fastq .

############
# run test
###########

# Set reference and parameters in nextflow.config
nextflow run nf_wochenende.nf  -with-timeline -with-report --fastq *R1.fastq 
