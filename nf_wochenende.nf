#!/usr/bin/env nextflow
/*
========================================================================================
                         nf_wochenende
========================================================================================
 Short and long read metagenomic alignment pipeline in Nextflow. Requires a fastq read file and a bwa indexed fasta reference genome 

 Colin Davenport, Lisa Hollstein

 #### Homepage / Documentation Changelog

v0.1.1  
v0.1.0  
v0.0.9  
v0.0.8  
v0.0.7  
v0.0.6  
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
params.metagenome = ""
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
    println "Version 0.0.5 by Colin Davenport and Lisa Hollstein with many further contributors"

    // File inputs
    //just R1 linked into dir
    //input_fastq = Channel.fromPath(params.fastq, checkIfExists: true)
    // both, but separate workflow dirs
    //input_fastq = Channel.fromPath("*_R{1,2}.fastq", checkIfExists: true)
    // Read inputs, SE read inputs should be possible
    input_fastq_R1 = Channel.fromPath("*_R1.fastq", checkIfExists: true)
    //input_fastq_R2 = Channel.fromPath("*_R2.fastq", checkIfExists: false)

    chunksize = Channel.value(1000)


    // run processes
   
    wochenende(input_fastq_R1)
    //wochenende(input_fastq_R1, input_fastq_R2)

    // run reporting
    reporting(wochenende.out.bam_txts.flatten())

    // run haybaler
    haybaler(reporting.out.us_csvs.collect())

    // create heattrees from haybaler output
    // needs R server
    //heattrees(haybaler.out.haybaler_heattree_csvs)

    // create heatmaps from haybaler ouput
    // needs R server
    //heatmaps(haybaler.out.haybaler_csvs.flatten())

    // run plots on the calmd_bams only
    plots(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais)

    // run growth_rate prediction software
    // growth_rate(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais, wochenende.out.bam_txts)

    // generate alignment stats
    //bam_stats(wochenende.out)

    // convert bam to cram format
    //convert_bam_cram(sort_bam.out)
    
    // multiqc
    //multiqc(bam_stats.out.collect(), bam_stats.out)

    // metagen window filter
    // metagen_window(wochenende.out.all, plots.out.window_txt)


}




/*
 *  Run wochenende
 *  Parcels the python script into a single Nextflow process
 *  Output - sorted bams for each step, and bam.txt files with read counts per chromosome.
 */

process wochenende {

    cpus = 16
	// If job fails, try again with more memory
	//memory { 40.GB * task.attempt }
    memory 40.GB
	errorStrategy 'terminate'

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
    //file fastq2


    output:
    //file "${prefix}*s.bam"
    //file "${prefix}*s.bam.bai"
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
    path "*", emit: all
    //path "ref.tmp", emit: ref_tmp
    

    script:
    name = fastq
    //prefix = fastq.name.toString().tokenize('.').get(0)
    String[] array
    array = fastq.name.toString().split('_R1');
    prefix = array[0]
    fastq_R2 = prefix + "_R2.fastq"
    println "Derived FASTQ R2 from R1 as: " + fastq_R2
    println params.metagenome
    println params.WOCHENENDE_DIR

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


    //export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}
    //export HAYBALER_DIR=${params.HAYBALER_DIR}

    //cp ${params.WOCHENENDE_DIR}/get_wochenende.sh .
    //bash get_wochenende.sh     

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

    ln -s ${launchDir}/$fastq_R2 .
    python3 run_Wochenende.py --metagenome ${params.metagenome} --threads $task.cpus --aligner $params.aligner $params.abra $params.mq --remove_mismatching $params.mismatches --readType $params.readType $params.prinseq $params.no_duplicate_removal --force_restart $fastq

    """

}


/*
 * Run reporting
 */

process reporting {
    cpus = 16

    conda params.conda_wochenende

    publishDir path: "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    file bamtxt

    output:
    path "*csv", emit: csvs
    path "*.rep.us.csv", emit: us_csvs
    path "*.rep.s.csv", emit: s_csvs

    script:

    """
    export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}

    cp ${params.WOCHENENDE_DIR}/reporting/basic_reporting.py .

    python3 basic_reporting.py --input_file $bamtxt --reference /mnt/ngsnfs/seqres/metagenref/bwa/2021_12_human_bact_arch_fungi_vir.fa --sequencer illumina --output_name $bamtxt
    """
}


/*
 * Run Haybaler
 * Requires Haybaler to be installed
 */

process haybaler {

    cpus = 12

    conda params.conda_haybaler

    publishDir path: "${params.outdir}/reporting", mode: params.publish_dir_mode

    input:
    file us_csv

    output:
    path "*haybaler*.csv", emit: haybaler_csvs
    path "*haybaler.csv", emit: haybaler_heattree_csvs
    path "haybaler_output"

    script:

    """
    cp ${params.HAYBALER_DIR}/haybaler.py .
    cp ${params.HAYBALER_DIR}/csv_to_xlsx_converter.py .
    cp ${params.WOCHENENDE_DIR}/haybaler/run_haybaler.sh .

    bash run_haybaler.sh
    """
}


/*
 * Run Heattrees
 */

process heattrees {
    cpus = 12

    conda params.conda_haybaler

    publishDir path: "${params.outdir}/reporting/haybaler_output", mode: params.publish_dir_mode

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
    cpus = 12

    conda params.conda_haybaler

    publishDir path: "${params.outdir}/reporting/haybaler_output", mode: params.publish_dir_mode

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

    cpus = 2
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
        publishDir path: "${params.outdir}/plots", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith('images')) "$filename"
                          else filename
                    }
    }
    


    input:
    file bam
    file bai
    //file fastq
    //file bam_txt


    output:
    file "images"
    path "*.calmd_cov_window.txt", emit: window_txt
    

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    cp -R ${params.WOCHENENDE_DIR}/plots/ .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    cp scripts/*.sh .
    bash runbatch_sambamba_depth.sh
    
    echo "INFO: Completed Sambamba depth and filtering"

    echo "INFO: Started Wochenende plot"
    cd plots
    cp ../*_window.txt . 
    cp ../*_window.txt.filt.csv .
    bash runbatch_wochenende_plot.sh >/dev/null 2>&1
    cd $launchDir
    echo "INFO: Completed Wochenende plot"


    """
}



/*
 *  Run growth rate
 */

process growth_rate {

    cpus = 2
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
    file "fit_results"
    

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    // run growth_rate scripts from current directory to avoid linking and output problems
    """
    cp -R ${params.WOCHENENDE_DIR}/growth_rate/ .
    cp -R ${params.WOCHENENDE_DIR}/scripts/ .
    cp scripts/*.sh .


    echo "INFO: Started bacterial growth rate analysis"
    cd $launchDir
    cp growth_rate/* .
        
    bash runbatch_bed_to_csv.sh  >/dev/null 2>&1 
        
    bash run_reproduction_determiner.sh  >/dev/null 2>&1
     
    cd $launchDir
    echo "INFO: Completed bacterial growth rate analysis, see growth_rate/fit_results/output for results"


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

    cpus = 8
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
    cpus = 2
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
    cpus = 2

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
