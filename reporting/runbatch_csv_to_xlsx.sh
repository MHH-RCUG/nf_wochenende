#!/bin/bash
set -eo pipefail
shopt -s nullglob

# Fabian Charly Friedrich and Colin Davenport
# Start from a folder containing csv files 
# bash runbatch_csv_to_xlsx.sh

# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)


# wochenende output files
for file in ./*rep.*.csv
do
        echo "$file"
	      python3 csv_to_xlsx_converter.py "$file" &
	wait
done
