# nf-snvcalling: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [bcftools mpileup and call](#bcftools) - The first mpileup part generates genotype likelihoods at each genomic position with coverage. The second call part makes the actual calls.
- [Annotation](#annotation) - Several databases like gnomAD, dbSNP, mirBASE, COSMIC and ExAC as well as the sequencing regions like selfchains, enhangers, repeat regions or mappability beds are used to create annotations to the variants. Annovar tool is used to annotate and make classifications.
- [Artifact Filtering](#filtering) - If runArtifactFilter is on, creates bias files for alternative and reference base and read positions and qualities, and plots error and basescore bias.
- [Deep Annotation](#annotation) - Several databases like gnomAD, dbSNP, mirBASE, COSMIC and ExAC as well as the sequencing regions like selfchains, enhangers, repeat regions or mappability beds are used to create annotations to the variants. Annovar tool is used to annotate and make classifications.
- [Filtering](#filtering) - Filtering can be applied to the annotated files for no-control samples.
- [SNV Extraction](#filtering) - Functional and somatic SNVs are extracted.
- [Plots](#visualization) - Several plots including rainfall, mutation distance, context frequencies, error plots, basescore distrubutions are created.
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Variant Calling

<details markdown="1">
<summary>Output files</summary>

- `metaid/`
  - `*snv_metaid.vcf.gz`: Raw variants.

</details>

### Annotation

<details markdown="1">
<summary>Output files</summary>

- `metaid/`
  - `metaid.deepanno.vcf.gz`: Annotated variants.

</details>

### Filtration

<details markdown="1">
<summary>Output files</summary>

- `metaid/`
  - `metaid_somatic_functional_snvs_conf_*_to_10.vcf`: Functional somatic snvs filtered through the confidence score
  - `metaid_somatic_snvs_conf_*_to_10.vcf`: Somatic snvs filtered through the confidence score
  - `metaid_somatic_ncRNA_snvs_conf_*_to_10.vcf`: Somatic ncRNA snvs filtered through the confidence score
  - `metaid_germline_functional_snvs_conf_*_to_10.vcf`: Functional germline snvs filtered through the confidence score
  - `metaid_reference_allele_base_qualities`: Allele base qualities of reference
  - `metaid_reference_allele_read_qualities`: Allele read qualities of reference
  - `metaid_alternative_allele_base_qualities`: Allele base qualities of variant file
  - `metaid_alternative_allele_read_qualities`: Allele read qualities of variant file
  - `metaid_base_score_bias_before_filter.pdf`: Base score error plot after bias filtration before filtration
  - `metaid_base_score_bias_after_filter_once.pdf`: First round of base score error plot after bias filtration
  - `metaid_sequence_specific_error_plot_before_filter.pdf`: Sequence spesific error plot after bias filtration before filtration
  - `metaid_sequence_specific_error_plot_after_filter_one.pdf`: First round of sequence spesific error plot after bias filtration
  - `metaid_sequencing_specific_error_plot_before_filter.pdf`: Sequencing spesific error plot after bias filtration before filtration
  - `metaid_sequencing_specific_error_plot_after_filter_one.pdf`: First round of sequencing spesific error plot after bias filtration

</details>

### Plots

<details markdown="1">
<summary>Output files</summary>

- `metaid/`
  - `metaid_allSNVdiagnosticPlots`: Merge of all diagnostic plots on quality of the variant.

</details>

### Reports

<details markdown="1">
<summary>Output files</summary>

- `metaid/`
  - `metaid_QC_values`: Statistics on QC.
  - `metaid_purityESTs`: Purity statistics.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
