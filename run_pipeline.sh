#!/bin/bash
module load  nextflow/22.07.1-edge
nextflow run main.nf -profile test_refgenie_37,singularity --outdir ref_37