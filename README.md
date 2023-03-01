
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<p align="center">
    <img title="nf-platypusindelcalling workflow" src="docs/images/nf-snvcalling.png" width=70%>
</p>
<p align="right">
    <img title="GHGA" src="docs/images/GHGA_short_Logo_orange.png" width=20%>
</p>
<p align="right">
    <img title="denbi" src="docs/images/denbi.png" width=20%>
</p>



## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-snvcalling** is a bioinformatics best-practice analysis pipeline from ODCF-OTP SNV Calling pipeline for somatic samples. It calls SNVs from both germline and somatic samples using bcftools mpileup, compares and filter outs germline spesific ones with samtools mpileup compare. This workflow uses various annotations from different databases, and applies broad filters accordingly.  Extensive QC plots serves functionality for high functional somatic mutation prioritization. 

For now, this workflow is only optimal to work in ODCF Cluster. The config file (conf/dkfz_cluster.config) can be used as an example. Running Annotation, DeepAnnotation and Filter steps are optinal and can be turned off using [runsnvAnnotation, runSNVDeepAnnotation, runSNVVCFFilter] parameters sequentialy.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 

**Important Notice:** The whole workflow is only ready for DKFZ cluster users for now, It is strongly recommended to them to read whole documentation before usage. This workflow works better with nextflow/22.07.1-edge in the cluster, It is recommended to use >22.07.1. Only SNV Calling part can be used for outside users, reference file and chromsome length file must be given for this.

## Pipeline summary

The pipeline has 5 main steps: SNV calling using mpileup, basic annotations, deep annotations, and filtering. Annotation and filtering steps are embeded with many plot generations. 


1. SNC Calling: 
   
   Bcftools mpileup ([`Bcftools mpileup`](https://samtools.github.io/bcftools/bcftools.html))
   : Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files. This is based on the original samtools mpileup command (with the -v or -g options) producing genotype likelihoods in VCF or BCF format, but not the textual pileup output.

2. Basic Annotations (--runSNVAnnotation True):

   In-house scripts to annotate with several databases like gnomAD and dbSNP.

   ANNOVAR ([`Annovar`](https://annovar.openbioinformatics.org/en/latest/))
   : annotate_variation.pl is used to annotate variants. The tool makes classifications for intergenic, intogenic, nonsynoymous SNP, frameshift deletion or large-scale duplication regions.
   
   Reliability and confidation annotations: It is an optional ste for mapability, hiseq, selfchain and repeat regions checks for reliability and confidence of those scores.

   Sequence ad Sequencing based error plots: Provides insights on predicted somatic SNVs.  

3. Deep Annotation (--runSNVDeepAnnotation True): 
   
   If basic annotations are applied, an extra optional step for number of extra indel annotations like enhancer, cosmic, mirBASE, encode databases can be applied too.

4. Filtering and Visualization (--runSNVVCFFilter True): 

   It is an optional step. Filtering is only required for the tumor samples with no-control and filtering can only be applied if basic annotation is performed. 

   SNV Extraction and Visualizations: SNVs can be extracted by certain minimum confidence level

   Visualization and json reports: Extracted SNVs are visualized and analytics of SNV categories are reported as JSON.

5. MultiQC (--skipmultiqc False):

   Produces pipeline level analytics and reports. 


**Please read** [usage](https://github.com/kubranarci/nf-snvcalling/blob/main/docs/usage.md)  before you start your won analysis. 

 
## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/))

3. Download [Annovar](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and set-up suitable annotation table directory to perform annotation. Example: 

 ```console
annotate_variation.pl -downdb wgEncodeGencodeBasicV19 humandb/ -build hg19
 ```

4. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   git clone https://github.com/kubranarci/nf-snvcalling.git
    ```

  before run do this to bin directory, make it runnable!:

  ```console
  chmod +x bin/*
  ```

   ```console
   nextflow run main.nf -profile test,YOURPROFILE --outdir <OUTDIR> --input <SAMPLESHEET>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker` and `singularity` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

5. Simple test run

   ```console
   nextflow run main.nf --outdir results -profile singularity,test
   ``` 

6. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run main.nf --input samplesheet.csv --outdir <OUTDIR> -profile <docker/singularity> --config test/institute.config
   ```
   
## Samplesheet columns

**sample**: The sample name will be tagged to the job

**tumor**: The path to the tumor file

**control**: The path to the control file, if there is no control will be kept blank.

## Data Requirements

All VCF and BED files need to be indexed with tabix and should be in the same folder!

[This section is for further]


## Documentation

The nf-psnvcalling pipeline comes with documentation about the pipeline [usage](https://github.com/kubranarci/nf-snvcalling/blob/main/docs/usage.md) and [output](https://github.com/kubranarci/nf-snvcalling/blob/main/docs/output.md).

## Credits

nf-snvcalling was originally written by Kuebra Narci kuebra.narci@dkfz-heidelberg.de.

The pipeline is originally written in workflow management language Roddy. [Inspired github page](https://github.com/DKFZ-ODCF/SNVCallingWorkflow)

We thank the following people for their extensive assistance in the development of this pipeline:

**TODO**
<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/platypusindelcalling for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
