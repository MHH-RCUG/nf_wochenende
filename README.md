
# nf_wochenende - a nextflow implementation of the Wochenende pipeline

**See documentation on our Wiki** at https://github.com/MHH-RCUG/nf_wochenende/wiki

This is a portable version of the alignment pipeline Wochenende which uses Nextflow. This should allow other users to start Wochenende more easily and efficiently.

* Colin Davenport
* Lisa Hollstein

# Status - mostly complete August 2022. 
* Graphical steps in R are least well tested and might need to be rerun manually on a configured R server


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


