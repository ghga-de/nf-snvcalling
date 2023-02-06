set -euo pipefail

usage() { echo "Usage: $0 [-p prefix] [-i vcf_file] [-t tumor] [-v basequal] [-b background] [-o output] [-r force] [-m median_threshold] [-w title] [-sp skipplots] [-rb refBaseQual] [-ab altBaseQual] [-ar altReadPos] [-rr refReadPos]" 1>&2; exit 1; }

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
			vcf=$2
			shift # past argument
	    	shift # past value
			;;
		-t)
			FILENAME_TUMOR_BAM=$2
			shift # past argument
	    	shift # past value
			;;
		-v)
			basequal=$2
			shift # past argument
	    	shift # past value
			;;
		-b)
			background=$2
			shift # past argument
	    	shift # past value
			;;
		-r)
			force=$2
			shift # past argument
	    	shift # past value
			;;
		-o)
			output=$2
			shift # past argument
	    	shift # past value
			;;
		-m)
			median_threshold=$2
			shift # past argument
	    	shift # past value
			;;
		-w)
			title=$2
			shift # past argument
	    	shift # past value
			;;

		-sp)
			skipplots=$2
			shift # past argument
	    	shift # past value
			;;
		-rb)
			refBaseQual=$2
			shift # past argument
	    	shift # past value
			;;
		-ab)
			altBaseQual=$2
			shift # past argument
	    	shift # past value
			;;
		-ar)
			altReadPos=$2
			shift # past argument
	    shift # past value
			;;
		-rr)
			refReadPos=$2
			shift # past argument
	    shift # past value
			;;
	esac
done

SEQUENCE_CONTEXT_COLUMN_INDEX=`cat ${vcf} | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "SEQUENCE_CONTEXT"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`

cat ${vcf} | perl -ne 'chomp; my $line=$_; if (/DP4=(\d+),(\d+),(\d+),(\d+);/) {my $fR=$1; my $rR=$2; my $fA=$3; my $rA=$4; my $MAF=($fA+$rA)/($fR+$rR+$fA+$rA); print "$line\t$MAF\n";} else { if (/^#CHROM/) { print "$line\tMAF\n";} else {print "$line\n";} };' >${prefix}.somatic_vcf.withMAF.vcf

MAF_COLUMN_INDEX=`cat ${prefix}.somatic_vcf.withMAF.vcf| grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, $_); my $columnIndex = first { $colnames[$_] eq "MAF"} 0..$#colnames; $columnIndex += 1; print "$columnIndex\n";'`

tripletBased_BQDistribution_plotter.R -v ${prefix}.somatic_vcf.withMAF.vcf -m '' -a ${FILENAME_TUMOR_BAM} \
        -p ${prefix} -b ${background} -o ${output} -R ${force} -c 1 \
        -f ${median_threshold} -s ${SEQUENCE_CONTEXT_COLUMN_INDEX} --MAFColumnIndex ${MAF_COLUMN_INDEX} -i 1 \
        -t "${title}" --skipPlots ${skipplots} \
        --refBaseQual ${refBaseQual} --altBaseQual ${altBaseQual} \
        --altReadPos ${altReadPos} --refReadPos ${refReadPos}
