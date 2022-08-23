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


############
# run nf_wochenende
###########

# Set reference and parameters in nextflow.config
nextflow run nf_wochenende.nf  -with-timeline -with-report --fastq *R1.fastq 
