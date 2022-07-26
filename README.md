
# nf_wochenende - a nextflow implementation of the Wochenende pipeline

**See documentation on our Wiki** at https://github.com/MHH-RCUG/nf_wochenende/wiki

This is a more portable version of the alignment pipeline Wochenende which uses Nextflow. This should allow other users to start Wochenende more easily and efficiently.

* Colin Davenport
* Lisa Hollstein

# Work in progress (Spring-Summer 2022). Do not use 


## Developer notes

Running nf_wochenende
* You **do not** need to edit and run setup.sh any more, this will be removed
* Edit config in `config.yaml`
* Edit config in `nextflow.config`  (scheduler, path to Wochenende, etc)
* Get a small read fastq file as input (present in the repo for testing)
* Edit `start_nf.sh`
* Run `bash start_nf.sh`

Output
* By default in output/wochenende

If missing or 0 sized bams in output
* Check most recent work folder
* cd work && ls -lrt
* Go into most recent folder and subfolder
* ls -l 
* cat .command.err
* It could be that trimming files or perldup have not be located correctly


Not touched to date
* shell scripts, wochenende_postprocess etc
* tidy all scripts to subdirs

Note that that only process used in the run_nf_wochenende.nf file is wochenende (see the workflow section). Other processes might be needed later, or may be removed.

![flowchart](https://user-images.githubusercontent.com/6094884/180654822-5d6d3129-e50a-485b-9fc8-52def87cb4b9.png)


