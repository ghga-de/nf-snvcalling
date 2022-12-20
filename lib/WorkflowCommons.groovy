//
// This file holds several functions common to the multiple workflows in the nf-core/viralrecon pipeline
//

class WorkflowCommons {

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "============================================================================="
            System.exit(1)
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Function to check whether primer BED file has the correct suffixes as provided to the pipeline
    //
    public static void checkPrimerSuffixes(primer_bed_file, primer_left_suffix, primer_right_suffix, log) {
        def total = 0
        def left  = 0
        def right = 0
        primer_bed_file.eachLine { line ->
            total += 1
            def name = line.split('\t')[3]
            if (name.contains(primer_left_suffix)) {
                left += 1
            } else if (name.contains(primer_right_suffix)) (
                right += 1
            )
        }
        if (total != (left + right)) {
            log.warn "=============================================================================\n" +
                "  Please check the name field (column 4) in the file supplied via --primer_bed.\n\n" +
                "  All of the values in that column do not end with those supplied by:\n" +
                "      --primer_left_suffix : $primer_left_suffix\n" +
                "      --primer_right_suffix: $primer_right_suffix\n\n" +
                "  This information is required to collapse the primer intervals into amplicons\n" +
                "  for the coverage plots generated by the pipeline.\n" +
                "==================================================================================="
        }
    }

    //
    // Function to get column entries from a file
    //
    public static ArrayList getColFromFile(input_file, col=0, uniqify=false, sep='\t') {
        def vals = []
        input_file.eachLine { line ->
            def val = line.split(sep)[col]
            if (uniqify) {
                if (!vals.contains(val)) {
                    vals << val
                }
            } else {
                vals << val
            }
        }
        return vals
    }

    //
    // Function that returns the number of lines in a file
    //
    public static Integer getNumLinesInFile(input_file) {
        def num_lines = 0
        input_file.eachLine { line ->
            num_lines ++
        }
        return num_lines
    }

    //
    // Function to generate an error if contigs in BED file do not match those in reference genome
    //
    public static void checkContigsInBED(fai_contigs, bed_contigs, log) {
        def intersect = bed_contigs.intersect(fai_contigs)
        if (intersect.size() != bed_contigs.size()) {
            def diff = bed_contigs.minus(intersect).sort()
            log.error "=============================================================================\n" +
                "  Contigs in primer BED file do not match those in the reference genome:\n\n" +
                "  ${diff.join('\n  ')}\n\n" +
                "  Please check:\n" +
                "    - Primer BED file supplied with --primer_bed\n" +
                "    - Genome FASTA file supplied with --fasta\n" +
                "============================================================================="
            System.exit(1)
        }
    }

    //
    // Function to read in all fields into a Groovy Map from Nextclade CSV output file
    //
    // See: https://stackoverflow.com/a/67766919
    public static Map getNextcladeFieldMapFromCsv(nextclade_report) {
        def headers   = []
        def field_map = [:]
        nextclade_report.readLines().eachWithIndex { row, row_index ->
            def vals = row.split(';')
            if (row_index == 0) {
                headers = vals
            } else {
                def cells = headers.eachWithIndex { header, header_index ->
                    def val = (header_index <= vals.size()-1) ? vals[header_index] : ''
                    field_map[header] = val ?: 'NA'
                }
            }
        }
        return field_map
    }

    //
    // Function to get number of variants reported in BCFTools stats file
    //
    public static Integer getNumVariantsFromBCFToolsStats(bcftools_stats) {
        def num_vars = 0
        bcftools_stats.eachLine { line ->
            def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
            if (matcher) num_vars = matcher[0][1].toInteger()
        }
        return num_vars
    }
}