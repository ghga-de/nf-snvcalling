# ghga-de/nf-snvcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.5dev

### `Fixed`

## v2.0.4

### `Added`

- nf-schema for parameter and samplesheet validation, fix `fasta.contain` logic [#48](https://github.com/ghga-de/nf-snvcalling/pull/49)
- `contig_file` parameter (`.bed` or `.tsv`) to restrict the analysis to specific contigs
- `process_annotation` resource label (8 CPUs, 12 GB, 32 h) for ANNOTATE_VCF, ANNOTATION_PIPES, ANNOVAR and FILTER_PEOVERLAP; replaces the overrides in the DKFZ cluster configs

### `Fixed`

- SM tag extraction [#48](https://github.com/ghga-de/nf-snvcalling/pull/49)
- AD field error [#47](https://github.com/ghga-de/nf-snvcalling/pull/49)
- Error detection in Perl scripts [#43](https://github.com/ghga-de/nf-snvcalling/issues/43), [#55](https://github.com/ghga-de/nf-snvcalling/pull/55): read errors are checked after each read loop and failed `zcat`/`tabix` pipes are detected on close, so a truncated or unreadable input stops the script instead of producing partial output
- `annotate_vcf.pl`: the a-file is opened once (with `-a file.vcf` the header lines were processed as variants)
- `snv_extractor_v1.pl`: the header is found when `--region` is used
- `findConfColumn.pl`: fails if there is no `#CHROM` header or no `CONFIDENCE` column (before, it printed a column index past the end of the line)

### `Removed`

- Email options and templates from the nf-core template (including `max_multiqc_email_size`)
- Unused Perl scripts (not called by any module or script):
  - `bin/confidenceAnnotation_SNVs.pl`, `bin/confidenceAnnotation_SNVsNoGermline.pl`, `bin/confidenceREAnnotation_SNVs.pl` (replaced by `bin/confidenceAnnotation_SNVs.py`)
  - `bin/mpileupOutputSplitter.pl`, `bin/vcfAmbiguousIndelStripper.pl`

### `Changed`

- CI: workflows run on `main` (was `master`), Java 17 and Prettier 3.9.0 are pinned, linting fixes

## v2.0.2 - 12.11.2024

### `Fixed`

- add cram/crai files to input list (assets/schema_input.json)
- add only a warning about that indices do not exists in converttovcf.json (bin/convertToStdVCF.py)
- change default arg for ENSEMBL VEP (conf/modules.config)
- upgrade samtools version for cram handling (modules/local/grep_samplename.nf)
- edit nextflow version to latest (nextflow.config)

## v2.0.1 - 09.10.2024

### `Added`

- nf-prov plugin is added.

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

- 20 output proper vcf by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/33
- re-arrange resources for dkfz cluster by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/34
- 35 raw vcf is not sorted by @kubranarci in https://github.com/ghga-de/nf-snvcalling/pull/36

**Full Changelog**: https://github.com/ghga-de/nf-snvcalling/compare/v1.0...v2.0.0

## v1.0dev - [date]

Initial release of ghga-de/nf-snvcalling, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`
