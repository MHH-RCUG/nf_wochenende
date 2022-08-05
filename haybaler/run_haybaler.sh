#!/bin/bash
# Lisa Hollstein, July 2022
# Based on the original script run_haybaler.sh by Sophia Poertner
# Run haybaler https://github.com/MHH-RCUG/haybaler/

version="0.23, August 2022"

#Args
readcount_limit=$1
rpmm_limit=$2

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
if [[ $count == 1 ]]
    then
      echo "ERROR: Haybaler needs 2 or more input CSV files to work. Use 2+ fastq files for the nf_wochenende pipeline or deselect haybaler in nf_wochenende.nf."
      exit 1
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

# actually run haybaler
python3 haybaler.py -i "$input_files" -p . -op $outputDir  -o haybaler.csv --readcount_limit $readcount_limit --rpmm_limit $rpmm_limit


# Move log file into log directory
mkdir -p $outputDir/logs
mv $outputDir/excluded_taxa.csv $outputDir/logs/
