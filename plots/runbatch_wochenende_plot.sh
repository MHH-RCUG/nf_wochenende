#!/bin/bash
# Colin Davenport, Jan 2020 - Nov 2021
# Run Wochenende plotting as batch, add SLURM


# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)

# Save output log in directory containing bams and preprocess script
output_log="plot_"$(date +%s)".log"
echo "INFO: output_log: " $output_log

# move existing images to a backup directory
images_backup="images_old_"$(date +%s)
mv images $images_backup

for i in $(ls *cov_window.txt.filt.csv)
        do
		# plot the prepared filtered csv files
		python wochenende_plot.py "$i" >> $output_log &

done
wait

