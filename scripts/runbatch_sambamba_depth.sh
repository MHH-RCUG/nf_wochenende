#!/bin/bash
# Erik Wessels, Colin Davenport Jan 2020 - Nov 2021
# Check window coverage on Wochenende sorted dup.bam output
# Use output for Python script to check coverage distribution


# Actually run for each BAM file
for i in `ls *calmd.bam`; do
	input=$i
	sec_input=${input%%.bam}
	#sec_in_bam=${input%%.bam}

	window=100000
	overlap=50000
	covMax=999999999
	
	# Get coverage depth in windows
	# x threads, Windows 100000, overlap 50000, -c minimum coverage. 
	#$path_sambamba depth window -t 8 --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &
	sambamba depth window -t 8 --max-coverage=$covMax --window-size=$window --overlap $overlap -c 0.00001 ${sec_input}.bam > ${sec_input}_cov_window.txt &

done
wait
