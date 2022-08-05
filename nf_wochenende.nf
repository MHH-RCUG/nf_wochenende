#!/usr/bin/env nextflow
/*
========================================================================================
                         nf_wochenende
========================================================================================
 Short and long read metagenomic alignment pipeline in Nextflow. Requires a fastq read file and a bwa indexed fasta reference genome 

 Colin Davenport, Lisa Hollstein

 #### Homepage / Documentation Changelog

v0.1.1  
v0.1.0  Raspir done, heat trees and heatmaps need to be manually tested as no R server in cluster
v0.0.9  Raspir integration underway
v0.0.8  Growth_rate fixed, plotting colours improved
v0.0.7  Plot (CD) and reporting (mainly Lisa) now fixed. Reporting fails if no data aligned to ref, fair enough.
v0.0.6  Add new mock reference seq and fast5 mock files for testing
v0.0.5  Remove get_wochenende.sh script functionality, still need WOCHENENDE_DIR defined in nextflow script for run_Wochenende.py, reporting semi-working, start metagen window filter
v0.0.4  First plot semi-working, start growth rate, test with bigger data
v0.0.3  Organize env variables, remove cluster submission bash code as now handled by nextflow
v0.0.2  Setup args
v0.0.1  init




----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is:
      conda activate nextflow
      nextflow run run_nf_wochenende.nf  --fasta /path/to/x.fa --fastq /path/x.fastq


    Arguments - all fully defined in script start.sh:
      --fasta [file]                  Path to Fasta reference. (Default: false)
      --fastq.gz [str]                fastq.gz read set
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, singularity

    Other
      --outdir [file]                 The output directory where the results will be saved (Default: './output')
    """.stripIndent()
}

// use modern nextflow
nextflow.enable.dsl = 2

/*
* Parameter defaults
*/

params.help=false
params.save_align_intermeds=true
params.outdir = "output"
params.publish_dir_mode = "copy"
params.fastq = ""
params.ref = ""
params.aligner = ""
params.mismatches = ""
params.nextera = ""
// params.abra = ""
params.mapping_quality = ""
params.readType = ""
params.debug = ""
params.longread = ""
params.no_dup_removal = ""
params.no_prinseq = ""
params.no_fastqc = ""
params.fastp = ""
params.trim_galore = ""




// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
* Workflow
*/


workflow {

    println "Starting run_nf_wochenende.nf"
    println "Version 0.1.0 by Colin Davenport, Tobias Scheithauer and Lisa Hollstein with many further contributors"

    // File inputs
    // R1 Read inputs, R2 reads are linked in by the process if they exist.
    input_fastq_R1 = Channel.fromPath("*_R1.fastq", checkIfExists: true)

    chunksize = Channel.value(1000)

    println "Using reference sequence: " + params.ref
    println "Using this WOCHENENDE_DIR: " + params.WOCHENENDE_DIR

    // Parameters - throw warnings at present
    if (params.mapping_quality != "") {
       params.mq = "--" + params.mapping_quality
    } else {
       params.mq = ""
    }

    if (params.no_abra) {
       params.abra = "--no_abra"
    } else {
       params.abra = ""
    }

    if (params.no_dup_removal) {
       params.no_duplicate_removal = "--no_duplicate_removal"
    } else {
       params.no_duplicate_removal = ""
    } 

    if (params.no_prinseq) {
       params.prinseq = "--no_prinseq"
    } else {
       params.prinseq = ""
    }


    // run processes

    // run main Wochenende process
    wochenende(input_fastq_R1)
    
    // run reporting
    reporting(wochenende.out.calmd_bam_txts.flatten())

    // run haybaler
    //haybaler(reporting.out.us_csvs.collect().flatten())
    haybaler(reporting.out.us_csvs.collect())

    // create heattrees from haybaler output
    // needs R server configured in config.yml
    //heattrees(haybaler.out.haybaler_heattree_csvs)

    // create heatmaps from haybaler ouput
    // needs R server
    //heatmaps(haybaler.out.haybaler_csvs.flatten())

    // run plots on the calmd_bams only
    plots(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais)

    // run growth_rate prediction step
    growth_rate(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais, wochenende.out.bam_txts)

    // run raspir steps
    raspir_fileprep(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais)

    //raspir(raspir_fileprep.out.collect())
    raspir(raspir_fileprep.out)

    // generate alignment stats
    //bam_stats(wochenende.out)

    // convert bam to cram format
    //convert_bam_cram(sort_bam.out)
    
    // multiqc
    //multiqc(bam_stats.out.collect(), bam_stats.out)

    // metagen window filter - @Lisa - I think this step is covered in plots and or growth rate so might not be needed?
    // metagen_window(wochenende.out.all, plots.out.window_txt)


}




/*
 *  Run wochenende
 *  Parcels the run_Wochenende.py python script into a single Nextflow process
 *  Output - sorted bams for each step, and bam.txt files with read counts per chromosome.
 *  Terminates on error, since this step provides data for all further steps.
 */

process wochenende {

    cpus = 16
	// If job fails, try again with more memory
	memory { 40.GB * task.attempt }
    //memory 40.GB
	//errorStrategy 'terminate'
    errorStrategy 'retry'

    // Use conda env defined in nextflow.config file
    // TODO - make a singularity container
    //conda '/home/hpc/davenpor/programs/miniconda3/envs/wochenende/'
    conda params.conda_wochenende


    tag "$name"
    label 'process_medium'
    
       
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/wochenende", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.bam')) "$filename"
                          else if (filename.endsWith('.bai')) "$filename"
                          else if (filename.endsWith('.bam.txt')) "$filename"
                          else if (filename.endsWith('.txt')) "$filename"
                          else if (filename.endsWith('.fastq')) "$filename"
                          else filename
                    }
    }
    


    input:
    file fastq


    output:
    path "*.bam", emit: bams
    path "*.s.bam", emit: s_bams
    path "*.calmd.bam", emit: calmd_bams
    path "*.mm.bam", emit: mm_bams
    path "*.dup.bam", emit: dup_bams
    path "*.bam.bai", emit: bam_bais
    path "*.s.bam.bai", emit: s_bam_bais
    path "*.calmd.bam.bai", emit: calmd_bam_bais
    path "*.mm.bam.bai", emit: mm_bam_bais
    path "*.dup.bam.bai", emit: dup_bam_bais
    path "*.bai"
    path "*.fastq"
    path "*.bam.txt", emit: bam_txts
    path "*.calmd.bam.txt", emit: calmd_bam_txts
    path "*.*", emit: all
    

    script:
    name = fastq
    //prefix = fastq.name.toString().tokenize('.').get(0)
    String[] array
    array = fastq.name.toString().split('_R1');
    prefix = array[0]
    fastq_R2 = prefix + "_R2.fastq"
    if (params.readType == "PE") {
        println "Derived FASTQ R2 from R1 as: " + fastq_R2
    }



    """
    export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}
    export HAYBALER_DIR=${params.HAYBALER_DIR}
    cp -f ${params.WOCHENENDE_DIR}/*.py .
    cp -f ${params.WOCHENENDE_DIR}/*.sh .
    cp -f ${params.WOCHENENDE_DIR}/*.config .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    cp -R ${params.WOCHENENDE_DIR}/reporting/ .
    cp -R ${params.WOCHENENDE_DIR}/dependencies/*.pl .
    cp scripts/*.sh .

    if [[ $params.readType == "PE" ]]
        then
        echo "Trying to link in R2, the second pair of the paired end reads"
        ln -s ${launchDir}/$fastq_R2 .
    fi
    python3 run_Wochenende.py --ref ${params.ref} --threads $task.cpus --aligner $params.aligner $params.abra $params.mq --remove_mismatching $params.mismatches --readType $params.readType $params.prinseq $params.no_duplicate_removal --force_restart $fastq

    """

}


/*
 * Run reporting
 */

process reporting {
    cpus = 1

    conda params.conda_wochenende
    errorStrategy 'ignore'
    //errorStrategy 'terminate'
    
    tag "$name"
    label 'process_medium'
    publishDir path: "${params.outdir}/reporting", mode: params.publish_dir_mode
	

    input:
    file bamtxt

    output:
    path "*csv", emit: csvs
    path "*.rep.us.csv", emit: us_csvs
    path "*.rep.s.csv", emit: s_csvs

    script:
    name = bamtxt

    """
    export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}

    cp ${params.WOCHENENDE_DIR}/reporting/basic_reporting.py .

    python3 basic_reporting.py --input_file $bamtxt --reference ${params.ref} --sequencer illumina --output_name $bamtxt
    """
}


/*
 * Run Haybaler
 * Requires Haybaler to be installed
 */

process haybaler {

    cpus = 1

    conda params.conda_haybaler
	//errorStrategy 'ignore'
    errorStrategy 'terminate'
    
    tag "$name"
    label 'process_medium'

    publishDir path: "${params.outdir}/haybaler", mode: params.publish_dir_mode

    input:
    file us_csv

    output:
    path "haybaler_output/*haybaler*.csv", emit: haybaler_csvs
    path "haybaler_output/*haybaler.csv", emit: haybaler_heattree_csvs
    path "haybaler_output/logs"

    script:
    name = "haybaler_input"

    // full run haybaler moved here to allow easy parameter changes
    // # Only run for *bam*.csv if files exist in current dir
    // # Only run for *bam*.txt if files exist in current dir
    // # Use  --readcount_limit 1 --rpmm_limit 10 for pipeline testing, use nextflow.config

    """/bin/bash
    cp ${params.HAYBALER_DIR}/haybaler.py .
    cp ${params.HAYBALER_DIR}/csv_to_xlsx_converter.py .
    cp ${params.WOCHENENDE_DIR}/haybaler/run_haybaler.sh .

    bash run_haybaler.sh ${params.haybaler_readcount_limit} ${params.haybaler_rpmm_limit}




    """
}


/*
 * Run Heattrees
 */

process heattrees {
    cpus = 1

    conda params.conda_haybaler
	errorStrategy 'ignore'
    //errorStrategy 'terminate'

    publishDir path: "${params.outdir}/haybaler", mode: params.publish_dir_mode

    input:
    file heattree_files

    output:
    path 'heattree_plots'
    path '*.csv'

    script:

    """
    cp ${params.WOCHENENDE_DIR}/haybaler/run_haybaler_tax.sh .
    cp ${params.HAYBALER_DIR}/haybaler_taxonomy.py .

    bash run_haybaler_tax.sh

    cp ${params.WOCHENENDE_DIR}/haybaler/run_heattrees.sh .
    cp ${params.HAYBALER_DIR}/create_heattrees.R .

    bash run_heattrees.sh
    """
}


/*
 * Run Heatmaps
 */

process heatmaps {
    cpus = 1

    conda params.conda_haybaler
	errorStrategy 'ignore'
    //errorStrategy 'terminate'

    publishDir path: "${params.outdir}/haybaler", mode: params.publish_dir_mode

    input:
    file heatmap_file

    output:
    path 'top*taxa/*'
    path '*filt.heatmap.csv'

    script:

    """
    cp ${params.WOCHENENDE_DIR}/runbatch_heatmaps.sh .
    cp ${params.HAYBALER_DIR}/create_heatmap.R .

    bash runbatch_heatmaps.sh
    """
}


/*
 *  Run plots
 */

process plots {

    cpus = 1
	// If job fails, try again with more memory if retry set
	memory { 8.GB * task.attempt }
	errorStrategy 'terminate'
    //errorStrategy 'ignore'
    //errorStrategy 'retry'

    // Use conda env defined in nextflow.config file
    conda params.conda_wochenende

    tag "$name"
    label 'process_medium'
    
       
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/plots", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('R1')) "$filename"
                          else filename
                    }
    }
    


    input:
    file bam
    file bai


    output:
    path "plots/images/*"
    path "*.calmd_cov_window.txt", emit: window_txt
    

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    cp -R ${params.WOCHENENDE_DIR}/plots/ .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    cp scripts/*.sh .
    bash runbatch_sambamba_depth.sh
    bash runbatch_metagen_window_filter.sh
    echo "INFO: Completed Sambamba depth and filtering"

    echo "INFO: Started Wochenende plot"
    cd plots
    cp ../*_window.txt . 
    cp ../*_window.txt.filt.csv .
    bash runbatch_wochenende_plot.sh >/dev/null 2>&1
    
        
    echo "INFO: Completed Wochenende plot"


    """
}



/*
 *  Run growth rate
 */

process growth_rate {

    cpus = 1
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	//errorStrategy 'terminate'
    errorStrategy 'ignore'
    //errorStrategy 'retry'

    // Use conda env defined in nextflow.config file
    conda params.conda_wochenende

    tag "$name"
    label 'process_medium'
    
       
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/growth_rate", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('fit_results')) "$filename"
                          else filename
                    }
    }
    


    input:
    file bam
    file bai
    file bam_txt


    output:
    //file "growth_rate/*"
    file "fit_results/*"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    // run growth_rate scripts from current directory to avoid linking and output problems
    """
    cp -R ${params.WOCHENENDE_DIR}/growth_rate/ .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    cp scripts/*.sh .


    echo "INFO: Started bacterial growth rate analysis"
    cp growth_rate/* .
        
    bash runbatch_bed_to_csv.sh  >/dev/null 2>&1 
        
    bash run_reproduction_determiner.sh  >/dev/null 2>&1
     
    echo "INFO: Completed bacterial growth rate analysis, see growth_rate/fit_results/output for results"


    """
}




/*
 *  Run raspir file preparation 
 */

process raspir_fileprep {

    cpus = 8
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	//errorStrategy 'terminate'
    errorStrategy 'ignore'
    //errorStrategy 'retry'

    // Use conda env defined in nextflow.config file
    conda params.conda_haybaler

    tag "$name"
    label 'process_medium'
    


    input:
    file bam
    file bai

    output:
    path "*.raspir.csv"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    cp -R ${params.WOCHENENDE_DIR}/raspir/ .

    echo "INFO: Started raspir analysis"
    cp raspir/* .

    bash run_SLURM_file_prep.sh $bam >/dev/null 2>&1
         
    echo "INFO: Completed raspir module"

    """
  
}


/*
 *  Run raspir
 */

process raspir {

    cpus = 1
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	//errorStrategy 'terminate'
    errorStrategy 'ignore'
    //errorStrategy 'retry'

    // Use conda env defined in nextflow.config file
    conda params.conda_haybaler

    tag "$name"
    label 'process_medium'
    
       
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/raspir", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('*.csv')) "$filename"
                          else filename
                    }
    }
    


    input:
    file input_csv
    //each file input_csv

    output:
    path "*.csv"

    script:
    prefix = input_csv.name.toString().tokenize('.').get(0)
    name = input_csv

    """
    cp -R ${params.WOCHENENDE_DIR}/raspir/ .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    
    echo "INFO: Started raspir analysis"
    cp raspir/* .

    python raspir.py $input_csv ${prefix}.csv >/dev/null 2>&1
    echo "INFO: Completed raspir"

    """
  
}



/*
 *  Convert BAM to coordinate sorted BAM, make stats, flagstat, idxstats
 */

process sort_bam {

    cpus = 8
	// If job fails, try again with more memory
	memory { 32.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.flagstat')) "$filename" 
                          else if (filename.endsWith('.idxstats')) "$filename" 
                          else if (filename.endsWith('.stats')) "$filename" 
                          else filename
                    }
    }


    input:
    file bam


    output:
    file "${prefix}.sorted.bam"
    file "${prefix}.sorted.bam.bai"
    file "${prefix}.sorted.bam.flagstat"
    file "${prefix}.sorted.bam.idxstats"
    file "${prefix}.sorted.bam.stats"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $name $bam
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}



/*
 *  make bam stats, flagstat, idxstats
 */

process bam_stats {

    cpus = 1
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.flagstat')) "$filename" 
                          else if (filename.endsWith('.idxstats')) "$filename" 
                          else if (filename.endsWith('.stats')) "$filename" 
                          else filename
                    }
    }


    input:
    file bam
    file bai
    file fastq
    file bam_txt



    output:
    file "${prefix}.sorted.bam.flagstat"
    file "${prefix}.sorted.bam.idxstats"
    file "${prefix}.sorted.bam.stats"

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
    samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
    samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
    """
}


/*
 *  Convert BAM to coordinate sorted BAM, make stats, flagstat, idxstats
 */

process convert_bam_cram {

    cpus = 8
	// If job fails, try again with more memory
	memory { 32.GB * task.attempt }
	errorStrategy 'retry'

    conda '/home/hpc/davenpor/programs/miniconda3/envs/wochenende/'


    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/samtools", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.cram')) "$filename"                         
                          else filename
                    }
    }


    input:
    file bam
    file bai
    file flagstat
    file idxstats
    file stats

    output:
    file "${prefix}.cram"
    file "${prefix}.cram.crai"


    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    // cram conversion  samtools view -@ {threads} -C -T {params.reference} -o {output.cram} {output.bam}

    """
    samtools view -@ $task.cpus -C -T $params.fasta -o ${prefix}.cram $bam
    samtools index ${prefix}.cram
    """
}






/*
 *  Multiqc
 */

process multiqc {
    cpus = 1
	// If job fails, try again with more memory
	memory { 4.GB * task.attempt }
	errorStrategy 'retry'

    // TODO - singularity 
    conda '/home/hpc/davenpor/programs/miniconda3/envs/bioinf/'
    
    tag "$name"
    label 'process_medium'
    if (params.save_align_intermeds) {
        publishDir path: "${params.outdir}/multiqc", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('.html')) "$filename"
                          else filename
                    }
    }

    input:
    path multiqc_files
    file flagstat
    file idxstats
    file stats

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    //when:
    //task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    name = "All stats files"
    
    """
    multiqc -f .
    
    """

}


process metagen_window {
    cpus = 1

    conda params.conda_wochenende

    publishDir path: "${params.outdir}", mode: params.publish_dir_mode

    input:
    file all
    file window_txt

    output:
    path "*window.txt.filt.csv", emit: window_filt
    path "*window.txt.filt.sort.csv", emit: window_sort

    script:

    """
    cp ${params.WOCHENENDE_DIR}/scripts/runbatch_metagen_window_filter.sh .
    bash runbatch_metagen_window_filter.sh
    """

}
