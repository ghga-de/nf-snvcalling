#!/bin/bash
module load  nextflow/22.10.7
nextflow run main.nf -profile dkfz_cluster_hg38,singularity --outdir /omics/odcf/analysis/evaluation_projects/workflow-validation/KN_sequencing/7_dec_sv_validations/WGS/ --input assets/samplesheet_hg38_WGS.csv --runcontigs "ALT_HLA" -resume