
set -euo pipefail

usage() { echo "Usage: $0 [-i vcf_file] [-p prefix] [-t mismatch_threshold] [-gnmax crit_gnomad_genomes_maxmaf] [-gemax crit_gnomad_exomes_maxmaf] [-lcmax crit_localcontrol_maxmaf] [-kgmax crit_1kgenomes_maxmaf] [-c isNoControlWorkflow] [-tb tumorbamfullpath] [-mq mapqual] [-r reference] [-o filenameSNVVCF] [-bq basequal]" 1>&2; exit 1; }

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
			mismatch_threshold=$2
			shift # past argument
	    shift # past value
			;;
		-gnmax)
			crit_gnomad_genomes_maxmaf=$2
			shift # past argument
	    shift # past value
			;;
		-gemax)
			crit_gnomad_exomes_maxmaf=$2
			shift # past argument
	    shift # past value
			;;
		-lcmax)
			crit_localcontrol_maxmaf=$2
			shift # past argument
	    shift # past value
			;;
		-kgmax)
			crit_1kgenomes_maxmaf=$2
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
		-mq)
			mapqual=$2
			shift # past argument
	    shift # past value
			;;
		-r)
			REFERENCE_GENOME=$2
			shift # past argument
	    shift # past value
			;;
		-o)
			filenameSNVVCF=$2
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

declare -r filenameSNVVCFTemp=${filenameSNVVCF}.tmp
declare -r filenameReferenceAlleleBaseScores=${filenameSNVVCF}_reference_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleBaseScores=${filenameSNVVCF}_alternative_allele_base_qualities.txt.gz
declare -r filenameAlternativeAlleleReadPositions=${filenameSNVVCF}_alternative_allele_read_positions.txt.gz
declare -r filenameReferenceAlleleReadPositions=${filenameSNVVCF}_reference_allele_read_positions.txt.gz

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

if [[ ${ref_spec-true} == "true" ]]
then
    hg38_params="--gnomAD_WGS_maxMAF=${crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${crit_localcontrol_maxmaf} --1000genome_maxMAF=${crit_1kgenomes_maxmaf}"
else
    hg38_params=""
fi


zcat ${vcf} | filter_PEoverlap.py ${noControlFlag} --alignmentFile=${tumorbamfullpath}\
	--mapq=$mapqual \
	--baseq=$basequal \
	--qualityScore=phred \
	--maxNumberOfMismatchesInRead=${mismatch_threshold--1} \ 
	--altBaseQualFile=${filenameAlternativeAlleleBaseScores}_NP \
	--refBaseQualFile=${filenameReferenceAlleleBaseScores}_NP \
	--altBasePositionsFile=${filenameAlternativeAlleleReadPositions}_NP \
	--refBasePositionsFile=${filenameReferenceAlleleReadPositions}_NP \
	--referenceFile=${REFERENCE_GENOME} \
		| confidenceAnnotation_SNVs.py ${noControlFlag} -i - -t 500 \
        	${hg38_params} > ${filenameSNVVCFTemp}

exitCode=$?
[[ $exitCode == 0 ]] && [[ -f ${filenameSNVVCFTemp} ]] && mv ${filenameSNVVCFTemp} ${filenameSNVVCF}
[[ $exitCode != 0 ]] && echo "SNV confidenceAnnotation with germline pipe returned non-zero exit code; temp file ${filenameSNVVCFTemp} not moved back" && exit 21

wait ${zipAlternativeAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleBaseScores" && exit 31
wait ${zipReferenceAlleleBaseScores} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleBaseScores" && exit 32
wait ${zipAlternativeAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipAlternativeAlleleReadPositions" && exit 33
wait ${zipReferenceAlleleReadPositions} ; [[ $? -gt 0 ]] && echo "Error from zipReferenceAlleleReadPositions" && exit 34

rm ${filenameAlternativeAlleleBaseScores}_NP ${filenameReferenceAlleleBaseScores}_NP ${filenameAlternativeAlleleReadPositions}_NP ${filenameReferenceAlleleReadPositions}_NP


bgzip -f ${filenameSNVVCF} && tabix -f -p vcf ${filenameSNVVCF}.gz
[[ $? != 0 ]] && echo "Error in creation of bgzipped vcf file and tabix index for it" && exit 41
