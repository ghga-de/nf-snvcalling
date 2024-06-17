#!/bin/bash
module load nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_hg38,singularity --outdir results2 --input assets/samplesheet_hg38_WGS.csv -resume