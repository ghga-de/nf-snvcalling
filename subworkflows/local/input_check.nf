//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK (samplesheet)
        .csv
        .splitCsv ( header:true, sep:',' )
        .map{ create_bam_channel(it) }
        .set {ch_sample}

    emit:
    ch_sample // channel: [ val(meta), [tumor,tumor.bai],[ control, control.bai]]
    versions = SAMPLESHEET_CHECK.out.versions
}

// Function to get list of [ sample, [ tumor],[ tumor_index], [control ], [control_index] ]
def create_bam_channel(LinkedHashMap row) {
// create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.iscontrol    = row.iscontrol

    // add path(s) of the fastq file(s) to the meta map
    def bam_meta = []

        if (!file(row.tumor).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Tumor file does not exist!\n${row.tumor}"
        }
        if (row.iscontrol) {
            if (!file(row.control).exists()) {
                if (row.control == 'dummy.bam') {
                    bam_meta = [  meta, file(row.tumor), file(row.tumor_index), [],[]  ]
                    meta.tumor_bam = file(row.tumor)
                    meta.tumor_bai = file(row.tumor_index)
                    meta.control_bam = []
                    meta.control_bai = []

                }
                else {
                    exit 1, "ERROR: Please check input samplesheet -> Control file does not exist!\n${row.control}"
                    }
            }

            bam_meta = [ meta, file(row.tumor), file(row.tumor_index), file(row.control), file(row.control_index) ]
            meta.tumor_bam = file(row.tumor)
            meta.tumor_bai = file(row.tumor_index)
            meta.control_bam = file(row.control)
            meta.control_bai = file(row.control_index)

        } else {
            bam_meta = [  meta, file(row.tumor), file(row.tumor_index), [],[]  ]
            meta.tumor_bam = file(row.tumor)
            meta.tumor_bai = file(row.tumor_index)
            meta.control_bam = []
            meta.control_bai = []
        }
    return bam_meta
}


