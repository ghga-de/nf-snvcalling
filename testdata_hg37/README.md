README

V_H021-1K9WSG sample from clinical data used to create test data.

samtools view --threads 4 -b /omics/odcf/analysis/evaluation_projects/workflow-validation/sequencing_full_test/whole_genome_sequencing/results_per_pid/V_H021-1K9WSG/alignment/tumor_V_H021-1K9WSG_merged.mdup.bam -L test.bed > tumor.bam 
samtools view --threads 4 -b /omics/odcf/analysis/evaluation_projects/workflow-validation/sequencing_full_test/whole_genome_sequencing/results_per_pid/V_H021-1K9WSG/alignment/buffy_coat_V_H021-1K9WSG_merged.mdup.bam -L test.bed > control.bam
samtools index -b tumor.bam
samtools index -b control.bam