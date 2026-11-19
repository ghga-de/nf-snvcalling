set -euo pipefail

usage() { echo "Usage: $0 [-p prefix] [-i filenameSomaticSnvs]  [-s filenameSomaticSnvsIndbSNP] [-t mincov] [-v minconfidencescore] [-b thatreshold] [-r rerun] [-j json_report]" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
  		-p)
			prefix=$2
			shift # past argument
	    	shift # past value
			;;
		-i)
			filenameSomaticSnvs=$2
			shift # past argument
	    	shift # past value
			;;
		-s)
			filenameSomaticSnvsIndbSNP=$2
			shift # past argument
	    	shift # past value
			;;
		-t)
			MINCOV=$2
			shift # past argument
	    	shift # past value
			;;
		-v)
			MIN_CONFIDENCE_SCORE=$2
			shift # past argument
	    	shift # past value
			;;
		-b)
			THA_SCORE_THRESHOLD=$2
			shift # past argument
	    	shift # past value
			;;
		-r)
			rerun=$2
			shift # past argument
	    	shift # past value
			;;
		-j)
			json_report=$2
			shift
			shift
			;;
	esac
done

# 2. MAF plots
# get mutant allele frequency ("MAF") for extracted somatic high confidence SNVs
makeMAFinput.pl ${filenameSomaticSnvs} "$MINCOV" "$MIN_CONFIDENCE_SCORE" > ${prefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt

[[ "$?" != 0 ]] && echo "There was a non-zero exit code in making the MAF input file" && exit 6

# count the obtained SNVs to output their number in the plot: < 50 will not be reliable!
snvnum=`grep -v "^#" ${filenameSomaticSnvs} | wc -l`
echo "snvnum is ${snvnum}"
snvindbSNP=` awk '{FS="\t"}{if(NR==2)print $5}'	${filenameSomaticSnvsIndbSNP}`
echo "snvindbSNP is ${snvindbSNP}"

# QC value $SNV_IN_DBSNP_RATIO will be written to $filenameQCvalues

export TEMP=$(mktemp -d)
export TMPDIR=$TEMP

# The MAF plot is made whenever there are SNVs, also when none of them is in dbSNP (#62):
# the plot then shows "0 in dbSNP" instead of an empty page.
# snvindbSNP is "NA" when the input has no DBSNP column (no dbSNP annotation).
if [ "$snvnum" != "0" ]; then
	MAF_plots.r ${prefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt ${snvnum} ${prefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf ${prefix} ${snvindbSNP}
	echo "MAF plots done"
	# QC value: fraction of SNVs in dbSNP
	if [[ "$snvindbSNP" =~ ^[0-9]+$ ]]; then
		SNV_IN_DBSNP_RATIO=$(awk -v a="$snvindbSNP" -v b="$snvnum" 'BEGIN { printf "%.5f", a / b }')
	else
		SNV_IN_DBSNP_RATIO="NA"
	fi
	echo "SNV_IN_DBSNP_RATIO is ${SNV_IN_DBSNP_RATIO}"
else
	SNV_IN_DBSNP_RATIO="NA"
	gs -sDEVICE=pdfwrite -o empty.pdf -c showpage
fi
rm -rf $TEMP
echo "check 1"

# infer baseQuality bias (PV4)-related THA score (QC value)
THA_SCORE=`determine_THA_score.R -i ${filenameSomaticSnvs}`
echo $THA_SCORE
[[ "$?" != 0 ]] && echo "There was a non-zero exit code in THA score determination script." && exit 24
# write the THA flag only for a numeric score above the threshold
# (before, [[ $(... | bc) ]] was always true, so the file was written for every sample)
if [[ "${THA_SCORE}" =~ ^[0-9.eE+-]+$ ]] && (( $(echo "${THA_SCORE} > ${THA_SCORE_THRESHOLD}" | bc -l) )); then
	echo -e "THA score\t${THA_SCORE}\n" >${prefix}_is_THA_affected.txt
fi

echo "check 2"

if [ "$json_report" = true ]; then
	# determine fraction of SNVs called as "synonymous SNV" among all exonic SNVs (QC value)
	EXONIC_CLASSIFICATION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "EXONIC_CLASSIFICATION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
	export EXONIC_CLASSIFICATION_COLUMN_INDEX=$((${EXONIC_CLASSIFICATION_COLUMN_INDEX}-1))
	ANNOVAR_FUNCTION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "ANNOVAR_FUNCTION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
	export ANNOVAR_FUNCTION_COLUMN_INDEX=$((${ANNOVAR_FUNCTION_COLUMN_INDEX}-1))
	SYNONYMOUS_RATIO=`grep -v '^#' ${filenameSomaticSnvs} | perl -F'\t' -ae 'BEGIN { my $total=0; my $synonymous=0; } if ($F[$ENV{"ANNOVAR_FUNCTION_COLUMN_INDEX"}] eq "exonic" ) {$total++; if ($F[$ENV{"EXONIC_CLASSIFICATION_COLUMN_INDEX"}] eq "synonymous SNV") {$synonymous++;}} END { print $total ? $synonymous/$total : "NA"; }'`
	echo "JSON REPORT WILL BE WRITTEN"
	# JSON number, or null for NA/empty values (#62)
	jnum() { if [[ "${1:-}" =~ ^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$ ]]; then echo "$1"; else echo null; fi; }
	echo -e "{" >${prefix}_QC_values${rerun}.json
	echo -e "\t\"snvnum\": $(jnum "${snvnum:-}")," >>${prefix}_QC_values${rerun}.json
	echo -e "\t\"snvindbSNP\": $(jnum "${snvindbSNP:-}")," >>${prefix}_QC_values${rerun}.json
	echo -e "\t\"snvInDbsnpRatio\": $(jnum "${SNV_IN_DBSNP_RATIO:-}")," >>${prefix}_QC_values${rerun}.json
	echo -e "\t\"synonymousRatio\": $(jnum "${SYNONYMOUS_RATIO:-}")," >>${prefix}_QC_values${rerun}.json
	echo -e "\t\"thaScore\": $(jnum "${THA_SCORE:-}")" >>${prefix}_QC_values${rerun}.json
	echo -e "}" >>${prefix}_QC_values${rerun}.json
	echo "done"
else
	echo "JSON REPORT WONT BE WRITTEN"
fi