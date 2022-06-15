
# nf_wochenende - a nextflow implementation of the Wochenende pipeline

**See documentation on our Wiki** at https://github.com/MHH-RCUG/nf_wochenende/wiki

Work in progress (Spring-Summer 2022).

We encountered portability problems when submitting jobs on non-SLURM clusters such as lsf etc. We therefore offer a nextflow wrapper to start Wochenende more easily and efficiently.

* Colin Davenport
* Lisa Hollstein


Running nf_wochenende
* Edit config in `config.yaml`
* Edit config in `nextflow.config`  (scheduler, etc)
* Get a small read fastq file as input
* Edit `start_nf.sh`
* Run `bash start_nf.sh`




Not touched to date
* shell scripts, postprocess etc

Note that that only process used in the run_nf_wochenende.nf file it wochenende (see the workflow section). Otherwise might be needed later, or may be removed.
