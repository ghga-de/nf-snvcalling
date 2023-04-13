#!/bin/bash
module load  nextflow/22.07.1-edge
nextflow run main.nf -profile dkfz_cluster_hg38,singularity --outdir results_38 --input seq2_testdata_snv/samplesheet_38.csv