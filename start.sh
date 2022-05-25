#!/bin/bash
# bash start.sh input.fastq

#ref=""
fastq=$1
#aligner="bwamem"
aligner="minimap2"

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
                -resume


done




