#!/usr/bin/env nextflow
/*
========================================================================================
                         nf_wochenende
========================================================================================
 Short and long read metagenomic alignment pipeline in Nextflow. Requires a fastq read file and a bwa indexed fasta reference genome 

 Colin Davenport 

 #### Homepage / Documentation Changelog

v0.1.1  
v0.1.0  
v0.0.9  
v0.0.8  
v0.0.7  
v0.0.6  
v0.0.5  
v0.0.4  
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
// params.no_dup_removal = ""
// params.no_prinseq = ""
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
    println "Version 0.0.3 by Colin Davenport and Lisa Hollstein with many further contributors"

    // File inputs
    //just R1 linked into dir
    //input_fastq = Channel.fromPath(params.fastq, checkIfExists: true)
    // both, but separate workflow dirs
    //input_fastq = Channel.fromPath("*_R{1,2}.fastq", checkIfExists: true)
    // Read inputs, SE read inputs should be possible
    input_fastq_R1 = Channel.fromPath("*_R1.fastq", checkIfExists: true)
    input_fastq_R2 = Channel.fromPath("*_R2.fastq", checkIfExists: false)

    chunksize = Channel.value(1000)


    // run processes
   
    
    wochenende(input_fastq_R1, input_fastq_R2)

    // run plots on the calmd_bams only
    // plots(wochenende.out.calmd_bams, wochenende.out.calmd_bam_bais)

    // generate alignment stats
    //bam_stats(wochenende.out)

    // convert bam to cram format
    //convert_bam_cram(sort_bam.out)
    
    // multiqc
    //multiqc(bam_stats.out.collect(), bam_stats.out)


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
    file fastq2


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
    path "*.bam.txt"
    

    script:
    //prefix = fastq.name.toString().tokenize('.').get(0)
    prefix = fastq.name.toString().tokenize('_').get(0)
    name = fastq
    //fastq_R2=prefix + "_R2.fastq"
    //print "Derived FASTQ R2 from R1 as: " + fastq_R2
    print params.metagenome
    print params.WOCHENENDE_DIR

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


    """
    export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}
    export HAYBALER_DIR=${params.HAYBALER_DIR}

    cp ${params.WOCHENENDE_DIR}/get_wochenende.sh .
    bash get_wochenende.sh        

    python3 run_Wochenende.py --metagenome ${params.metagenome} --threads $task.cpus --aligner $params.aligner $params.abra $params.mq --remove_mismatching $params.mismatches --readType $params.readType $params.prinseq $params.no_duplicate_removal --debug --force_restart $fastq

    """

}




/*
 *  Run plots
 */

process plots {

    cpus = 2
	// If job fails, try again with more memory
	memory { 8.GB * task.attempt }
	errorStrategy 'terminate'
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
    file "${prefix}.sam"
    

    script:
    prefix = bam.name.toString().tokenize('.').get(0)
    name = bam

    """
    export WOCHENENDE_DIR=${params.WOCHENENDE_DIR}
    export HAYBALER_DIR=${params.HAYBALER_DIR}

    cp ${params.WOCHENENDE_DIR}/get_wochenende.sh .
    bash get_wochenende.sh 

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
