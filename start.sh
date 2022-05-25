#!/bin/bash
# bash start.sh input.fastq

# Uncomment the following to use them in the nextflow submission script below
#ref=""
fastq=$1
#aligner="bwamem"
#aligner="minimap2short"
aligner="minimap2long"
#aligner="ngmlr"
#mismatches="--remove_mismatching 3"
mismatches="--remove_mismatching 250"
nextera=""
#nextera="--nextera"
#abra=""
abra="--no_abra"
#mapping_quality="--mq20"
mapping_quality="--mq30"
readType="SE"
#readType="PE"
debug="--debug"
longread="--longread"
no_dup_removal="--no_duplicate_removal"
no_prinseq="--no_prinseq"
no_fastqc="no_fastqc"
fastp="--fastp"
trim_galore="--trim_galore"


git pull -q 

# check args
if [ -z $fastq ]
        then
        echo ""
        echo "## Usage: Input fastq required. bash start.sh in.fastq ## "
        exit
fi

# Run script - Paired end reads R2 will be calculated by replacing R1 with R2
# Uncomment/adapt the only line you want to run

#2021_12, minor update for 2021_10 ref. 
ref="/mnt/ngsnfs/seqres/metagenref/bwa/2021_12_human_bact_arch_fungi_vir.fa"
#viruses
#ref="/mnt/ngsnfs/seqres/metagenref/bwa/EZV0_1_database2_cln.fasta"
#ref="/mnt/ngsnfs/seqres/metagenref/bwa/nci_viruses.fa"





echo "INFO: Arguments check OK, starting job submission "


for i in `ls *R1.fastq`

        do
                echo $i
                #sbatch run_Wochenende_SLURM.sh $i
                # run specifying fastq reads and fasta reference as arg
                nextflow run run_nf_wochenende.nf  -with-timeline -with-report -with-dag flowchart.dot \
                --fasta $ref \
                --fastq $fastq \
                --test yes  \
                --aligner $aligner \
                --mismatches $mismatches
                
                -resume


done

fastq=$1
#aligner="bwamem"
#aligner="minimap2short"
aligner="minimap2long"
#aligner="ngmlr"
#mismatches="--remove_mismatching 3"
mismatches="--remove_mismatching 250"
nextera=""
#nextera="--nextera"
#abra=""
abra="--no_abra"
#mapping_quality="--mq20"
mapping_quality="--mq30"
readType="SE"
#readType="PE"
debug="--debug"
longread="--longread"
no_dup_removal="--no_duplicate_removal"
no_prinseq="--no_prinseq"
no_fastqc="no_fastqc"
fastp="--fastp"
trim_galore="--trim_galore"






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
