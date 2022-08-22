
# nf_wochenende - a nextflow implementation of the Wochenende pipeline

**See documentation on our Wiki** at https://github.com/MHH-RCUG/nf_wochenende/wiki

This is a portable version of the metagenomic alignment pipeline Wochenende which uses Nextflow. This should allow most users to use Wochenende more easily and efficiently.


Wochenende runs alignment of short reads (eg Illumina) or long reads (eg Oxford Nanopore) against a reference sequence. It is relevant for genomics and metagenomics. Wochenende is simple (python script), portable and is easy to configure with a central config file. 



## Features

Features include (see programs listed below at the bottom of this page)
- QC (Fastqc)
- pre alignment duplicate removal (perldup)
- pre alignment poor sequence removal (Prinseq - used for single ended reads only)
- trimming (Trimmomatic or fastp or trim galore or ea-utils)
- alignment (bwa mem, or minimap2 for short reads, minimap2 or ngmlr for long reads)
- SAM-> BAM conversion (samtools and sambamba)
- AlignerBoost Mapping Quality recalculation (in testing April-June 2021)
- Report % aligned reads (samtools)
- Output unmapped reads as fastq (samtools)  (from v1.4) 
- Post-alignment duplicate removal (Samtools from v1.7.8, Sambamba)
- Removal reads with x mismatches (bamtools), adjustable from v1.7.3
- Realignment (Abra2)
- MD tag marking (Samtools)
- Normalization (to Bacteria per Human cell, RPMM Reads Per Million sequenced reads per Million reference bases etc, see Reporting below for details)
- Visualization (chromosome coverage, intended for bacteria in metagenomics projects) (from v1.4)
- Growth rate estimation. Estimate how fast detected bacteria are growing
- Rare species prediction (Raspir)

Project Haybaler https://github.com/MHH-RCUG/haybaler allows postprocessing of Wochenende results:
- collation/integration of multiple reports (reporting csv or bam.txt files) using Python Pandas
- prepare results for heatmaps
- create heatmaps using multiple different R libraries

### Did you know ? 
Wochenende means weekend in German. The original developer, Tobias, called the pipeline Wochenende, because you can start it running and go off to enjoy your weekend early (at least, that was the plan!).

## Preprint 

Please view and cite the Wochenende preprint at https://www.biorxiv.org/content/10.1101/2022.03.18.484377v2


## Installation

Quick install and run guide for nf_wochenende 
* See full installation docs here: https://github.com/MHH-RCUG/nf_wochenende/wiki
* Edit config in `config.yaml`
* Edit config in `nextflow.config`  (cluster scheduler or local server, path to Wochenende and haybaler git clones and conda envs)

* Get a small read fastq file as input (present in the repo for testing)
* Edit `start_nf.sh`
* Run `bash start_nf.sh`

Output
* By default in output/wochenende
* One subfolder is output for each successful stage (nextflow process)
* See docs here: https://github.com/MHH-RCUG/nf_wochenende/wiki




![flowchart](https://user-images.githubusercontent.com/6094884/180654822-5d6d3129-e50a-485b-9fc8-52def87cb4b9.png)


## Tools used

- [Alignerboost](https://github.com/Grice-Lab/AlignerBoost). GPL3, in dependencies folder.
- [ABRA2](https://github.com/mozack/abra2)
- [bamtools](https://github.com/pezmaster31/bamtools)
- [BWA](https://github.com/lh3/bwa)
- [ea-utils](https://github.com/ExpressionAnalysis/ea-utils)
- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FastUniq](https://sourceforge.net/projects/fastuniq/)
- [Minimap2](https://github.com/lh3/minimap2)
- [NGMLR](https://github.com/philres/ngmlr)
- [perldup](https://github.com/richardmleggett/scripts/blob/master/remove_pcr_duplicates.pl) Already copied to dependencies folder with permission. Developed by [Richard Leggett](https://github.com/richardmleggett).
- [PRINSEQ](http://prinseq.sourceforge.net/)
- [sambamba](https://github.com/biod/sambamba)
- [samtools](https://github.com/samtools/samtools)
- [trim_galore](https://github.com/FelixKrueger/TrimGalore)
- [trimmomatic](https://github.com/timflutre/trimmomatic)

Postprocessing
- [Haybaler](https://github.com/MHH-RCUG/haybaler)

Optional extras
- [growth rate](https://github.com/MHH-RCUG/Wochenende/tree/master/growth_rate)
- [raspir](https://github.com/mmpust/raspir)



## Nextflow version authors

* Colin Davenport
* Lisa Hollstein
* Ilona Rosenboom
