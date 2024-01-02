#!/bin/bash
module load  nextflow/22.10.7
nextflow run main.nf -profile dkfz_cluster_hg37,singularity --outdir /omics/odcf/analysis/evaluation_projects/workflow-validation/KN_sequencing/7_dec_sv_validations/ --input assets/samplesheet_val_WGS.csv -resume