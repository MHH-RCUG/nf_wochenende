#!/bin/bash
# Lisa Hollstein, July 2022
# Based on the original script run_haybaler.sh by Sophia Poertner
# Run haybaler https://github.com/MHH-RCUG/haybaler/

version="0.22, June 2021"

outputDir=haybaler_output
if [[ ! -d $outputDir ]]
then
    echo "INFO: Creating directory:" $outputDir
    mkdir $outputDir
fi


input_files=""

# Only run for *bam*.csv if files exist in current dir
count=$(ls -1 *.bam*.csv 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for csv in *.bam*.csv
    do
      input_files="$input_files;$csv"
    done
fi


# Only run for *bam*.txt if files exist in current dir
count=$(ls -1 *.bam*.txt 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for csv in *.bam*.txt
    do
      input_files="$input_files;$csv"
    done
fi

python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv
# for pipeline testing only!!
#python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit 1 --rpmm_limit 10

# Move log file into log directory
mkdir -p $outputDir/logs
mv $outputDir/excluded_taxa.csv $outputDir/logs/
