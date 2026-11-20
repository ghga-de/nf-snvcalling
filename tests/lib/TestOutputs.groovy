import java.security.MessageDigest
import java.util.zip.GZIPInputStream

/**
 * Helpers for tests/main.nf.test (loaded via libDir in nf-test.config).
 */
class TestOutputs {

    /** Samplesheet for the test BAMs, with absolute paths (nf-test runs from its own directory). */
    static String samplesheet(String outputDir, String baseDir) {
        def sheet = new File("${outputDir}/samplesheet.csv")
        sheet.parentFile.mkdirs()
        def bam = { String s -> "${baseDir}/testdata/WES_LL_${s}_1_test.bam" }
        sheet.text = "sample,tumor,tumor_index,control,control_index\n" +
                     "test,${bam('T')},${bam('T')}.bai,${bam('N')},${bam('N')}.bai\n"
        return sheet.absolutePath
    }

    /** Lines of a text file, gzipped or not. */
    static List<String> lines(File f) {
        def stream = f.name.endsWith('.gz') ? new GZIPInputStream(new FileInputStream(f)) : new FileInputStream(f)
        return stream.withReader { it.readLines() }
    }

    /**
     * md5 of the data lines only: VCF headers contain dates, versions and paths
     * that change between runs and machines, the variants do not.
     */
    static String bodyMd5(File f) {
        def md = MessageDigest.getInstance('MD5')
        lines(f).findAll { !it.startsWith('#') }.each { md.update((it + '\n').getBytes('UTF-8')) }
        return md.digest().encodeHex().toString()
    }

    /** Checks that do not depend on the annotation tool (regression checks for #45, #59, #60, #62). */
    static void checkOutputs(String outputDir) {
        def res = "${outputDir}/results/test"

        def somatic = new File("${res}/snvs_test_somatic_snvs_conf_8_to_10.vcf")
        assert somatic.exists()
        def somatic_lines = lines(somatic)
        def header = somatic_lines.find { it.startsWith('#CHROM') }.split('\t') as List
        def variants = somatic_lines.findAll { !it.startsWith('#') }
        assert variants.size() > 0

        // #45: every variant has a SEQUENCE_CONTEXT with 10 bases on each side
        def ctx = header.indexOf('SEQUENCE_CONTEXT')
        assert ctx >= 0
        def bad_ctx = variants.findAll { !(it.split('\t')[ctx] ==~ /[ACGTN]{10},[ACGTN]{10}/) }
        assert bad_ctx.isEmpty() : "malformed SEQUENCE_CONTEXT: ${bad_ctx.take(3)}"

        // standardised VCF
        def std_vcfs = new File("${res}/standard_vcf").listFiles()?.findAll { it.name.endsWith('.vcf.gz') } ?: []
        assert std_vcfs.size() > 0 : "no standardised VCF in ${res}/standard_vcf"
        std_vcfs.each { f ->
            def std_lines = lines(f)
            // #60: no whitespace in INFO values (VCF 4.2)
            def bad_info = std_lines.findAll { !it.startsWith('#') && it.split('\t')[7] =~ /\s/ }
            assert bad_info.isEmpty() : "whitespace in INFO of ${f.name}: ${bad_info.take(3)}"
            // #59: mappability is kept in INFO
            assert std_lines.any { it.startsWith('##INFO=<ID=MAP,') } : "no MAP INFO header in ${f.name}"
        }

        // #62: diagnostics plots are merged into one PDF
        assert new File("${res}/snvs_test_allSNVdiagnosticsPlots.pdf").exists()
    }

    /**
     * Stable outputs for the snapshot:
     *  - names of all published files of the sample (plots and tables are only listed: they contain dates/random layout)
     *  - md5 of the variant lines of every published VCF
     */
    static Map stableOutputs(String outputDir) {
        def res = new File("${outputDir}/results/test")
        def files = []
        res.eachFileRecurse(groovy.io.FileType.FILES) { files << res.toPath().relativize(it.toPath()).toString() }
        def vcf_md5 = [:]
        files.findAll { it.endsWith('.vcf') || it.endsWith('.vcf.gz') }.sort().each { vcf_md5[it] = bodyMd5(new File(res, it)) }
        return [ files: files.sort(), vcf_variant_md5: vcf_md5 ]
    }
}
