# ghga-de/nf-snvcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**Full Changelog**: https://github.com/ghga-de/nf-snvcalling/compare/v1.0...v2.0.0
## v1.0dev - [date]

Initial release of ghga-de/nf-snvcalling, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## v2.0.0 - 24.06.2024

### `Added`

- assets/config/convertToStdVCF.json and bin/convertToStdVCF.py 
    - Option to output VCF files (all) in standard format (4.2) is added. Also, TSV formatted confidence annotated/filtrated files are being converted into standard VCF.

- `MAFCommon` tag is added to the INFO column to mark the common/recurrent artefacts.

- Minor changes: 
    - output names of the VCF files.

- modules/local/sort_nonstandard_vcf.nf
    - sorted output for nonstandard raw vcf files. 
    
### `Fixed`

- Contig processing is only available for hg38 reference. ALT and/or HLA contigs can be given external in a file. 
    - Automatic generation of HLA/ALT contigs is now possible through tumor BAM instead of fasta.

- Conda links in nf-core modules are fixed. 
    - NOTE: Conda environments are not available for the pipeline. Holding conda environment.yml links the same in default creates an error even when enable_conda is false.  

- bin/confidenceAnnotation_SNVs.py
    - Flag parsing is generic now. 

- Better dealing with ALT and HLA contigs. 
   - modules/local/get_contigs.nf fixed.
   - ALT and HLA contig extraction is fixed 
   - Turn off error when there is no alignment to contigs.

- bin/vcf_pileup_compare_allin1_basecount.pl 
   - sorted results

- modules/local/seq_context_annotator.nf
   - sorted results
   
### `Dependencies`

### `Deprecated`

- FREQ-based filtering is removed from bin/confidenceAnnotation_SNVs.py.
- bcftools1.9 biocontainers docker changed with kubran/bcftools:v1.9

## What's Changed
* 20 output proper vcf by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/33
* re-arrange resources for dkfz cluster by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/34
* 35 raw vcf is not sorted by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/36


## v2.0.1 - 09.10.2024

### `Added`
 - nf-prov plugin is added.
