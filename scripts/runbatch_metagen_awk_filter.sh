#!/bin/bash
# Colin Davenport, April 2020 - Nov 2021
# Run multiqc report
# Collect mapping stats from flagstat
# Run filter: Keep all lines in the bam.txt where column 3 (reads aligned)
# is greater than X (here probably 20). Good for idxstats files i.e. bam.txt files from Wochenende


# Setup SLURM using data parsed from config.yaml
source $WOCHENENDE_DIR/scripts/parse_yaml.sh
eval $(parse_yaml $WOCHENENDE_DIR/config.yaml)

# Run samtools stats
echo "INFO:  Running samtools stats"
count=$(ls -1 *trm.s.bam 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
	for bam in `ls *trm.s.bam`
		do
		$path_samtools stats $bam > $bam.stats &
	done
fi
wait
# Count PE reads only if fix.s.bam files are present
count=$(ls -1 *fix.s.bam 2>/dev/null | wc -l)
if [[ $count != 0 ]]
    then
    for bam in `ls *fix.s.bam`
		do
		$path_samtools stats $bam > $bam.stats &
	done
fi
wait
# Count calmd bam files
count=$(ls -1 *calmd.bam 2>/dev/null | wc -l)
if [[ $count != 0 ]]
	then
	for bam in `ls *calmd.bam`
		do
		$path_samtools stats $bam > $bam.stats &
	done
fi
wait

# Run multiqc
if [[ "${USE_MULTIQC}" == "yes" ]]; then
	echo "INFO:  Running multiqc"
	multiqc -f .  2>&1 &
fi


# Collate mapping stats
out="mapped_percent.txt"
echo "INFO:  Generating Wochenende mapping stats to $out"
echo "Wochenende mapping stats" > $out
for file in `ls *flagstat.txt`
	do
	echo -n $file "\t"  >> $out
	grep "mapped (" $file >> $out
done


echo "INFO:  Filtering and sorting  Wochenende output bam.txt files"
## Filter: $3>= means column3 (no. reads assigned to taxon) must have 20 or more reads
# Sorted descending in column 3 for within experiment clarity 
for i in `find . -name "*.bam.txt"`
        do
		# Run directly, else can adversely impact slurmdbd database.
		#awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv &
		awk -F "\t" '$3>=20' $i | sort -t$'\t' -k3 -nr > $i.filt.sort.csv &
done
wait

echo "INFO:  Script runbatch_metagen_awk_filter.sh completed"
