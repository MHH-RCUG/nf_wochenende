#!/bin/bash
# Input data: requires tables output by haybaler_taxonomy.R https://github.com/MHH-RCUG/haybaler
# First run Wochenende, then wochenende_postprocess.sh. (or haybaler_taxonomy.sh manually)
# exclude GC, ref length, any host chr etc
# Sophia Poertner, Colin Davenport, 2020-2021
# Installation: First prepare server with R packages, see https://github.com/MHH-RCUG/haybaler
# Usage: bash run_heattrees.sh
# Adjustments for nf_wochenende by Lisa Hollstein

# Args
rscript_bin=$1

set_ulimits () {
        echo "INFO: trying to set ulimits higher from 8GB cstack to 32GB cstack"
        ulimit -s 32384
}

prepare_files () {
  echo "INFO: Preparing files for R heat-tree creation"
  for infile in {RPMM,bacteria_per_human_cell}*haybaler_taxa.csv
        do
        echo "Running on " "$infile"

        count_infile=$(ls -1 $infile 2>/dev/null | wc -l)
        if [[ $count_infile != 0 ]]
            then
            #exclude mouse, human, mito
            grep -v "^chr" "$infile" | grep -v "^1_1_1" > "${infile%_haybaler_taxa.csv}"_heattree_filt1.csv

            # using tab delimiters, cut a max of 2000 columns out excluding cols 2-3.
            cut -f1,4-2000 "${infile%_haybaler_taxa.csv}"_heattree_filt1.csv  > "${infile%_haybaler_taxa.csv}"_heattree.csv

            # cleanup
            rm "${infile%_haybaler_taxa.csv}"_heattree_filt1.csv
        fi
  done
}


create_heattrees () {
  rscript_bin="/usr/bin/Rscript"
  echo "INFO: Starting batch heat-tree creation"
    if [[ ! -f $rscript_bin ]]
            then
            echo "INFO: Rscript binary not found, aborting. Could not find this, is R installed? " $rscript_bin
            exit
    fi
    echo "INFO: Using rscript binary: " $rscript_bin


    # create heat-tree for each heatmap.csv file
    # check if correct files present, else do not run the script, which would create empty files.
    count=$(ls -1 {RPMM,bacteria_per_human_cell}*heattree.csv 2>/dev/null | wc -l)
    if [[ $count != 0 ]]
        then
        for heattreecsv in {RPMM,bacteria_per_human_cell}*heattree.csv
            do
            echo "INFO: Creating heat-tree for file: $heattreecsv"
            # run local
            $rscript_bin create_heattrees.R "$heattreecsv"
        done
    else
        echo "Could not find heat-tree input files of type heattree.csv in the directory"
    fi
}


# Actually run functions
set_ulimits
prepare_files
create_heattrees

echo "INFO: Cleanup - creating directories and moving files"
# check if directories exist, create them if not
mkdir -p heattree_plots
mkdir -p heattree_plots/RPMM_background_heattrees
mkdir -p heattree_plots/RPMM_no_background_heattrees
mkdir -p heattree_plots/bphc_background_heattrees
mkdir -p heattree_plots/bphc_no_background_heattrees

# move heat-trees in directories
RPMM_count_no_background_pdf=$(ls -1 RPMM_*no_background_heattree.pdf 2>/dev/null | wc -l)
if [[ $RPMM_count_no_background_pdf != 0 ]]
    then
    mv RPMM_*no_background_heattree.pdf heattree_plots/RPMM_no_background_heattrees
fi

RPMM_count_background_pdf=$(ls -1 RPMM_*background_heattree.pdf 2>/dev/null | wc -l)
if [[ $RPMM_count_background_pdf != 0 ]]
    then
    mv RPMM_*background_heattree.pdf heattree_plots/RPMM_background_heattrees
fi


bphc_count_no_background_pdf=$(ls -1 bacteria_per_human_cell*no_background_heattree.pdf 2>/dev/null | wc -l)
if [[ $bphc_count_no_background_pdf != 0 ]]
    then
    mv bacteria_per_human_cell*no_background_heattree.pdf heattree_plots/bphc_no_background_heattrees
fi

bphc_count_background_pdf=$(ls -1 bacteria_per_human_cell*background_heattree.pdf 2>/dev/null | wc -l)
if [[ $bphc_count_background_pdf != 0 ]]
    then
    mv bacteria_per_human_cell*background_heattree.pdf heattree_plots/bphc_background_heattrees
fi

count_all_samples=$(ls -1 *all_samples_heattree.pdf 2>/dev/null | wc -l)
if [[ $count_all_samples != 0 ]]
    then
    mv *all_samples_heattree.pdf heattree_plots
fi


# move log files
echo "moving empty_heattree_samples.txt to logs/"

mkdir -p logs

count_empty_heattree_sample=$(ls -1 empty_heattree_samples.txt 2>/dev/null | wc -l)
if [[ $count_empty_heattree_sample != 0 ]]
    then
    mv empty_heattree_samples.txt logs/
fi


echo "INFO: Heat tree script completed"

