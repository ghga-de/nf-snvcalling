set -euo pipefail

usage() { echo "Usage: $0 [-p prefix] [-i filenameSomaticSnvs]  [-s filenameSomaticSnvsIndbSNP] [-t mincov] [-v minconfidencescore] [-b thatreshold] [-r rerun]" 1>&2; exit 1; }

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

if [ "$snvindbSNP" != "0" ]; then
	if [ "$snvnum" != "0" ]; then
		MAF_plots.r ${prefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt ${snvnum} ${prefix}_MAF_conf_${MIN_CONFIDENCE_SCORE}_to_10.pdf ${prefix} ${snvindbSNP}
		[[ "$?" != 0 ]] && echo "There was a non-zero exit code in MAF plotting" && exit 7
		echo "MAF plots done"
	else
		gs -sDEVICE=pdfwrite -o empty.pdf -c showpage
	# calculate SNV_IN_DBSNP_RATIO snvindbSNP/snvnum (QC value)
	fi
	##SNV_IN_DBSNP_RATIO=$( (expr $snvindbSNP / $snvnum) | bc -l)
	SNV_IN_DBSNP_RATIO=$(echo "scale=5; $snvindbSNP / $snvnum" | bc)
	echo "SNV_IN_DBSNP_RATIO is ${SNV_IN_DBSNP_RATIO}"
else
	SNV_IN_DBSNP_RATIO="NA" # no output produced, don't include in "convert" later
	gs -sDEVICE=pdfwrite -o empty.pdf -c showpage
fi
rm -rf $TEMP
echo "check 1"

# infer baseQuality bias (PV4)-related THA score (QC value)
THA_SCORE=`determine_THA_score.R -i ${filenameSomaticSnvs}`
echo $THA_SCORE
[[ "$?" != 0 ]] && echo "There was a non-zero exit code in THA score determination script." && exit 24
[[ $(echo "${THA_SCORE} > ${THA_SCORE_THRESHOLD}" | bc -l) ]] && echo -e "THA score\t${THA_SCORE}\n" >${prefix}_is_THA_affected.txt

echo "check 2"
# determine fraction of SNVs called as "synonymous SNV" among all exonic SNVs (QC value)
EXONIC_CLASSIFICATION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "EXONIC_CLASSIFICATION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
export EXONIC_CLASSIFICATION_COLUMN_INDEX=$((${EXONIC_CLASSIFICATION_COLUMN_INDEX}-1))
ANNOVAR_FUNCTION_COLUMN_INDEX=`cat ${filenameSomaticSnvs} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "ANNOVAR_FUNCTION"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`
export ANNOVAR_FUNCTION_COLUMN_INDEX=$((${ANNOVAR_FUNCTION_COLUMN_INDEX}-1))
SYNONYMOUS_RATIO=`grep -v '^#' ${filenameSomaticSnvs} | perl -F'\t' -ae 'BEGIN { my $total=0; my $synonymous=0; } if ($F[$ENV{"ANNOVAR_FUNCTION_COLUMN_INDEX"}] eq "exonic" ) {$total++; if ($F[$ENV{"EXONIC_CLASSIFICATION_COLUMN_INDEX"}] eq "synonymous SNV") {$synonymous++;}} END { print $synonymous/$total; }'`
echo "check 3"
echo -e "{" >${prefix}_QC_values${rerun}.json
echo -e "\t\"snvnum\": ${snvnum:-NA}," >>${prefix}_QC_values${rerun}.json
echo -e "\t\"snvindbSNP\": ${snvindbSNP:-NA}," >>${prefix}_QC_values${rerun}.json
echo -e "\t\"snvInDbsnpRatio\": ${SNV_IN_DBSNP_RATIO:-NA}," >>${prefix}_QC_values${rerun}.json
echo -e "\t\"synonymousRatio\": ${SYNONYMOUS_RATIO:-NA}," >>${prefix}_QC_values${rerun}.json
echo -e "\t\"thaScore\": ${THA_SCORE:-NA}" >>${prefix}_QC_values${rerun}.json
echo -e "}" >>${prefix}_QC_values${rerun}.json
echo "done"
