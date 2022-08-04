#!/bin/bash
echo "bash script to start the nf_wochenende pipeline"
echo "Requires nextflow and bwa in path. Maybe do this first: conda activate wochenende"

# get most current version - for developers
git pull -q

# ignore system JAVA_HOME, use that supplied by wochenende conda env for nextflow
unset JAVA_HOME
unset WOCHENENDE_DIR
unset HAYBALER_DIR

# cleanup work/ 
#rm -rf work

# get test files, create reference (small files). Do it yourself for bigger test references 
# eg mock community from SRA https://github.com/colindaven/wochenende_manuscript/blob/main/mock/download_fastq.sh
# eg mock ref file human22_zymo_test.fa
#bwa index test/data/ref.fa
#cp -f test/data/*.fastq . && bwa index test/data/ref.fa
#cp -f test/data/*.fastq .

############
# run test
###########
#nextflow run nf_wochenende.nf  -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner bwamem --remove_mismatching 2 --mq30 --readType PE  --no_dup_removal --no_abra --fastq *R1.fastq

# lisa version with corrected args
#nextflow run nf_wochenende.nf  -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner bwamem --mismatches 2 --mapping_quality mq30 --readType PE  --no_prinseq --no_abra --fastq *R1.fastq

# PE with dup removal
#nextflow run nf_wochenende.nf  -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner bwamem --mismatches 2 --mapping_quality mq30 --readType PE --no_abra --fastq *R1.fastq

# Single ended reads, with --ref not --metagenome
#nextflow run nf_wochenende.nf  -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner bwamem --mismatches 2 --mapping_quality mq30 --readType SE  --no_abra --fastq *R1.fastq

# big reference
#nextflow run nf_wochenende.nf  -resume -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner bwamem --mismatches 2 --mapping_quality mq30 --readType SE  --no_abra --fastq *R1.fastq 

# big reference, long read
#nextflow run nf_wochenende.nf  -resume -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/seqres/metagenref/wochenende/2021_12_human_bact_arch_fungi_vir.fa --aligner minimap2long --mismatches 30 --mapping_quality mq30 --readType SE --longread --no_abra --fastq *R1.fastq 

#small ref, faster reporting
nextflow run nf_wochenende.nf  -resume -with-timeline -with-report --ref /mnt/beegfs/scratch/bioinformatics/colin/dev/nf_wochenende/current/nf_wochenende/test/data/ref.fa --aligner bwamem --mismatches 2 --mapping_quality mq30 --readType SE  --no_abra --fastq *R1.fastq 




#########
# Notes
###########

# run pipeline, no bash variables allowed as empty strings evaluate to true in groovy
#nextflow run nf_wochenende.nf  -with-timeline -with-report -with-dag flowchart.dot --metagenome Ath --aligner minimap2long --remove_mismatching 250 --mq30 --readType SE --longread --no_dup_removal --no_abra --fastq test_sm.fastq

# python3 /run_Wochenende.py --ref  --threads 16 --aligner minimap2-long   250 SE  --debug --force_restart test_sm.fastq



# All options
# --threads 12
#--longread (implies aligner is not bwa-mem, no prinseq, no duplicate removal)
#--aligner bwamem
#--aligner minimap2short   (for Illumina reads)
#--aligner minimap2long    (for nanopore reads)
#--aligner ngmlr
#--nextera  - remove Nextera adapters with Trimmomatic, not default Ultra II / Truseq adapters
#--no_abra  - no read realignment
#--mq20     - remove reads with a mapping quality of less than 20. Less stringent than MQ30, required for raspir https://github.com/mmpust/raspir
#--mq30     - remove reads with a mapping quality of less than 30
#--readType SE - single ended reads
#--readType PE - paired end reads
#--debug
#--force_restart
#--remove_mismatching 3 (remove those reads with 3 or more mismatches) # questionable for very long reads
#--remove_mismatching 250 (remove those reads with 250 or more mismatches) # for long reads set to ~10% of median read length.
#--no_duplicate_removal  - do not remove duplicate reads
#--no_prinseq   - do not filter out low complexity initial reads using prinseq (default in this file after 2020_11)
#--no_fastqc
#--testWochenende - runs the test scripts with test reads vs a testDB and checks if all seems well.
#--fastp - fastp is recommended as an alternative trimmer to Trimmomatic if you are having adapter problems
#--trim_galore - trim_galore is the best adapter trimmer for nextera reads
