
# nf_wochenende - a nextflow implementation of the Wochenende pipeline

**See documentation on our Wiki** at https://github.com/MHH-RCUG/nf_wochenende/wiki

This is a more portable version of the alignment pipeline Wochenende which uses Nextflow. This should allow other users to start Wochenende more easily and efficiently.

* Colin Davenport
* Lisa Hollstein

Work in progress (Spring-Summer 2022).


Running nf_wochenende
* Edit config in `config.yaml`
* Edit config in `nextflow.config`  (scheduler, etc)
* Get a small read fastq file as input
* Edit `start_nf.sh`
* Run `bash start_nf.sh`



Not touched to date
* shell scripts, wochenende_postprocess etc
* tidy all scripts to subdirs

Note that that only process used in the run_nf_wochenende.nf file is wochenende (see the workflow section). Other processes might be needed later, or may be removed.
