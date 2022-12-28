#!/bin/bash
module load  nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster,singularity --outdir result --input testdata_hg37/samplesheet_37.csv