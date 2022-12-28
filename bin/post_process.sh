
set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-p prefix] [-t runArtifactFilter] [-q generateExtendedQcPlots] [-c isNoControlWorkflow] [-tb tumorbamfullpath] [-ofpe filterPeOverlap_OPTIONS] [-ovb filterVcfForBias_OPTIONS] [-r reference] [-oc CONFIDENCE_OPTS] [-rs ref_spec] [-bq basequal]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			vcf=$2
			shift # past argument
	    shift # past value
			;;
		-p)
			prefix=$2
			shift # past argument
	    shift # past value
			;;
		-t)
			runArtifactFilter=$2
			shift # past argument
	    shift # past value
			;;
		-q)
			generateExtendedQcPlots=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			isNoControlWorkflow=$2
			shift # past argument
	    shift # past value
			;;
		-tb)
			tumorbamfullpath=$2
			shift # past argument
	    shift # past value
			;;
		-ofpe)
			filterPeOverlap_OPTIONS=$2
			shift # past argument
	    shift # past value
			;;
		-ovb)
			filterVcfForBias_OPTIONS=$2
			shift # past argument
	    shift # past value
			;;
		-r)
			REFERENCE_GENOME=$2
			shift # past argument
	    shift # past value
			;;
		-oc)
			CONFIDENCE_OPTS=$2
			shift # past argument
	    shift # past value
			;;
		-rs)
			ref_spec=$2
			shift # past argument
	    shift # past value
			;;
		-bq)
			basequal=$2
			shift # past argument
	    shift # past value
			;;
			
	esac
done

declare -r filenameSNVVCF=${prefix} # Cut off .gz from the end to get an intermediate file.
declare -r filenameSNVVCFTemp=${filenameSNVVCF}.tmp
declare -r filenamePCRerrorMatrix=${filenameSNVVCF}_sequence_error_matrix.txt
declare -r filenameSequencingErrorMatrix=${filenameSNVVCF}_sequencing_error_matrix.txt
declare -r filenameBiasMatrixSeqFile=${filenameSNVVCF}_sequence_specific_bias_matrix.txt
declare -r filenameBiasMatrixSeqingFile=${filenameSNVVCF}_sequencing_specific_bias_matrix.txt
declare -r filenameSomaticSNVsTmp=${filenameSNVVCF}_somatic_snvs_for_bias.vcf
declare -r filenameSequenceErrorPlotPreFilter=${filenameSNVVCF}_sequence_specific_error_plot_before_filter.pdf
declare -r filenameSequencingErrorPlotPreFilter=${filenameSNVVCF}_sequencing_specific_error_plot_before_filter.pdf
declare -r filenameBaseScoreBiasPreFilter=${filenameSNVVCF}_base_score_bias_before_filter.pdf
declare -r filenameSequenceErrorPlotTmp=${filenameSNVVCF}_sequence_specific_error_plot_after_filter_once.pdf
declare -r filenameSequencingErrorPlotTmp=${filenameSNVVCF}_sequencing_specific_error_plot_after_filter_once.pdf
declare -r filenameBaseScoreBiasTmp=${filenameSNVVCF}_base_score_bias_after_filter_once.pdf
declare -r filenameQCValues=${filenameSNVVCF}_QC_values.tsv
declare -r filenameReferenceAlleleBaseScores=${filenameSNVVCF}_reference_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleBaseScores=${filenameSNVVCF}_alternative_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleReadPositions=${filenameSNVVCF}_alternative_allele_read_positions.txt.gz
declare -r filenameReferenceAlleleReadPositions=${filenameSNVVCF}_reference_allele_read_positions.txt.gz


declare -r filenamePCRerrorMatrixFirst=${filenameSNVVCF}_sequence_error_matrix_first.txt
declare -r filenameSequencingErrorMatrixFirst=${filenameSNVVCF}_sequencing_error_matrix_first.txt
declare -r filenameBiasMatrixSeqFileFirst=${filenameSNVVCF}_sequence_specific_bias_matrix_first.txt
declare -r filenameBiasMatrixSeqingFileFirst=${filenameSNVVCF}_sequencing_specific_bias_matrix_first.txt

declare -r filenamePCRerrorMatrixSecond=${filenameSNVVCF}_sequence_error_matrix_second.txt
declare -r filenameSequencingErrorMatrixSecond=${filenameSNVVCF}_sequencing_error_matrix_second.txt
declare -r filenameBiasMatrixSeqFileSecond=${filenameSNVVCF}_sequence_specific_bias_matrix_second.txt
declare -r filenameBiasMatrixSeqingFileSecond=${filenameSNVVCF}_sequencing_specific_bias_matrix_second.txt

declare -r filenameControlMedian=${filenameSNVVCF}_control_median.txt

# create BaseScore FIFOs and their consumer processes (zip and write to target file)
# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}
mkfifo ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP
cat ${filenameAlternativeAlleleBaseScores}_NP | bgzip -f >${filenameAlternativeAlleleBaseScores} & zipAlternativeAlleleBaseScores=$!
cat ${filenameReferenceAlleleBaseScores}_NP | bgzip -f >${filenameReferenceAlleleBaseScores} & zipReferenceAlleleBaseScores=$!
cat ${filenameAlternativeAlleleReadPositions}_NP | bgzip -f >${filenameAlternativeAlleleReadPositions} & zipAlternativeAlleleReadPositions=$!
cat ${filenameReferenceAlleleReadPositions}_NP | bgzip -f >${filenameReferenceAlleleReadPositions} & zipReferenceAlleleReadPositions=$!

if [[ ${isNoControlWorkflow-false} == "false" ]]
then
    noControlFlag=""
else
    noControlFlag="--nocontrol"
fi
# If this is for the pancancer workflow, then also create a DKFZ specific file.
if [[ ${runArtifactFilter-true} == true ]]
then
	cat ${vcf} | filter_PEoverlap.py ${noControlFlag} --alignmentFile=${tumorbamfullpath} ${filterPeOverlap_OPTIONS} \
						--altBaseQualFile=${filenameAlternativeAlleleBaseScores}_NP --refBaseQualFile=${filenameReferenceAlleleBaseScores}_NP \
						--altBasePositionsFile=${filenameAlternativeAlleleReadPositions}_NP --refBasePositionsFile ${filenameReferenceAlleleReadPositions}_NP --referenceFile=${REFERENCE_GENOME} \
						| confidenceAnnotation_SNVs.py ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 0 -f ${filenameSomaticSNVsTmp} \
                            ${ref_spec} > ${filenameSNVVCFTemp}.tmp

    [[ $? != 0 ]] && echo "Error in first iteration of confidence annotation" && exit 2

	NRSOMSNV=`grep -v "^#" ${filenameSomaticSNVsTmp} | wc -l`
	echo -e "SOMATIC_SNVS_UNFILTERED\t${NRSOMSNV}">> ${filenameQCValues}

	mv ${filenameSNVVCFTemp}.tmp ${filenameSNVVCFTemp}

    wait ${zipAlternativeAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleBaseScores" && exit 31
    wait ${zipReferenceAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleBaseScores" && exit 32
    wait ${zipAlternativeAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleReadPositions" && exit 33
    wait ${zipReferenceAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleReadPositions" && exit 34
    rm ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP

	createErrorPlots.py --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequencingErrorPlotPreFilter} \
		--errorType=sequencing_specific --errorFile=${filenameSequencingErrorMatrix} --plot_title='Sequencing strand bias before guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in first creation of error matrix and plot (sequencing)" && exit 3

	createErrorPlots.py --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequenceErrorPlotPreFilter} \
		--errorType=sequence_specific --errorFile=${filenamePCRerrorMatrix} --plot_title='PCR strand bias before guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in first creation of error matrix and plot (sequence/PCR)" && exit 4

    if [[ ${generateExtendedQcPlots} == true ]]; then
        cp ${filenameSomaticSNVsTmp} ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering
        tripletBased_BQRatio_plotter.R -v ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering \
        	-r ${filenameReferenceAlleleBaseScores} -a ${filenameAlternativeAlleleBaseScores} -t ${basequal} -p "Differences" \
        	-o ${filenameBaseScoreBiasPreFilter} -d "Base Quality Bias Plot for PID ${prefix} before guanine oxidation filter" & plotBaseScoreBiasBeforeFiltering=$!
    fi

	filterVcfForBias.py --vcfFile=${filenameSNVVCFTemp} --referenceFile=${REFERENCE_GENOME} ${filterVcfForBias_OPTIONS} \
		--sequence_specificFile=${filenamePCRerrorMatrix} --sequencing_specificFile=${filenameSequencingErrorMatrix} \
		--bias_matrixSeqFile=${filenameBiasMatrixSeqFile} --bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile} --vcfFileFlagged=${prefix}_flagged.vcf \
			| confidenceAnnotation_SNVs.py ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 1 -f ${filenameSomaticSNVsTmp} \
        		${ref_spec} > ${filenameSNVVCFTemp}.tmp

	[[ $? != 0 ]] && echo "Error in first filtering and/or second interation of confidence annotation" && exit 5

	mv ${filenameSNVVCFTemp}.tmp ${filenameSNVVCFTemp}
	mv ${filenamePCRerrorMatrix} ${filenamePCRerrorMatrixFirst}
	mv ${filenameSequencingErrorMatrix} ${filenameSequencingErrorMatrixFirst}
	mv ${filenameBiasMatrixSeqFile} ${filenameBiasMatrixSeqFileFirst}
	mv ${filenameBiasMatrixSeqingFile} ${filenameBiasMatrixSeqingFileFirst}

	createErrorPlots.py --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequencingErrorPlotTmp} \
		--errorType=sequencing_specific --errorFile=${filenameSequencingErrorMatrix} --plot_title='Sequencing strand bias after first round of guanine oxidation filter'

    [[ $? != 0 ]] && echo "Error in second creation of error matrix and plot (sequencing)" && exit 6

	createErrorPlots.py --vcfFile=${filenameSomaticSNVsTmp} --referenceFile=NA --outputFile=${filenameSequenceErrorPlotTmp} \
		--errorType=sequence_specific --errorFile=${filenamePCRerrorMatrix} --plot_title='PCR strand bias after first round of guanine oxidation filter'

	[[ $? != 0 ]] && echo "Error in second creation of error matrix and plot (sequence/PCR)" && exit 7

    if [[ ${generateExtendedQcPlots} == true ]]; then
        cp ${filenameSomaticSNVsTmp} ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering
        tripletBased_BQRatio_plotter.R -v ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering -r ${filenameReferenceAlleleBaseScores} \
        	-a ${filenameAlternativeAlleleBaseScores} -t ${basequal} -p 'Differences' -o ${filenameBaseScoreBiasTmp} \
        	-d "Base Quality Bias Plot for PID ${PID} after first round of guanine oxidation filter" & plotBaseScoreBiasAfterFirstFiltering=$!
    fi

	filterVcfForBias.py --vcfFile=${filenameSNVVCFTemp}  --referenceFile=${REFERENCE_GENOME} --sequence_specificFile=${filenamePCRerrorMatrix} \
		--sequencing_specificFile=${filenameSequencingErrorMatrix} ${filterVcfForBias_OPTIONS} --bias_matrixSeqFile=${filenameBiasMatrixSeqFile} \
		--bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile} --vcfFileFlagged=${prefix}_flagged.vcf \
		| confidenceAnnotation_SNVs.py ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 2 \
        	${ref_spec} > ${filenameSNVVCF}

	[[ $? != 0 ]] && echo "Error in second filtering and/or third iteration of confidence annotation" && exit 8

	mv ${filenamePCRerrorMatrix} ${filenamePCRerrorMatrixSecond}
	mv ${filenameSequencingErrorMatrix} ${filenameSequencingErrorMatrixSecond}
	mv ${filenameBiasMatrixSeqFile} ${filenameBiasMatrixSeqFileSecond}
	mv ${filenameBiasMatrixSeqingFile} ${filenameBiasMatrixSeqingFileSecond}

	rm ${filenameSomaticSNVsTmp}
	rm ${filenameSNVVCFTemp}

	[[ $? != 0 ]] && echo "Error in moving the vcf file and index or in removing the temporary files" && exit 9

	if [[ ${generateExtendedQcPlots} == true ]]; then
    wait ${plotBaseScoreBiasBeforeFiltering} ; [[ $? -gt 0 ]] && echo "Error in first creation of base score bias plot" && exit 37
    [[ -f ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering ]] && rm ${filenameSomaticSNVsTmp}.forBSBPlotsBeforeFiltering
    wait ${plotBaseScoreBiasAfterFirstFiltering} ; [[ $? -gt 0 ]] && echo "Error in second creation of base score bias plot" && exit 38
    [[ -f ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering ]] && rm ${filenameSomaticSNVsTmp}.forBSBPlotsAfter1stFiltering
	fi
else
	cat ${vcf} | filterPeOverlap ${noControlFlag} --alignmentFile=${tumorbamfullpath} ${filterPeOverlap_OPTIONS} \
		--altBaseQualFile=${filenameAlternativeAlleleBaseScores}_NP --refBaseQualFile=${filenameReferenceAlleleBaseScores}_NP \
		--altBasePositionsFile=${filenameAlternativeAlleleReadPositions}_NP --refBasePositionsFile=${filenameReferenceAlleleReadPositions}_NP --referenceFile=${REFERENCE_GENOME} \
			| confidenceAnnotation_SNVs.py ${noControlFlag} -i - ${CONFIDENCE_OPTS} \
        	${ref_spec} > ${filenameSNVVCFTemp}

    exitCode=$?
    [[ $exitCode == 0 ]] && [[ -f ${filenameSNVVCFTemp} ]] && mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
    [[ $exitCode != 0 ]] && echo "SNV confidenceAnnotation with germline pipe returned non-zero exit code; temp file ${filenameSNVVCFTemp} not moved back" && exit 21

    wait ${zipAlternativeAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleBaseScores" && exit 31
    wait ${zipReferenceAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleBaseScores" && exit 32
    wait ${zipAlternativeAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleReadPositions" && exit 33
    wait ${zipReferenceAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleReadPositions" && exit 34

    rm ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP
fi

bgzip -f ${filenameSNVVCF} && tabix -f -p vcf ${filenameSNVVCF}.gz
[[ $? != 0 ]] && echo "Error in creation of bgzipped vcf file and tabix index for it" && exit 41
