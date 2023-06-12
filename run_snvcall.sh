#!/bin/bash
module load  nextflow/22.10.7
nextflow run main.nf -profile dkfz_cluster_hg38,singularity --outdir /omics/odcf/analysis/evaluation_projects/workflow-validation/KN_sequencing/10_may_sv_validations/ --input assets/samplesheet.csv -resume