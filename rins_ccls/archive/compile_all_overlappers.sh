#!/usr/bin/env bash
set -x


mapq=20
index="hg19"
core="bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned"
basedir=`pwd`

#	If passed 1 fast[aq], check for chimeric reads.
#	If passed 2 fast[aq], also check for anchors with paired read run.

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--mapq 20] [--index hg19]" 
	echo "[--core bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bowtie2.herv_k113.unaligned]"
	echo
	echo "Defaults:"
	echo "  mapq  ..... : $mapq"
	echo "  index ..... : $index"
	echo "  core  ..... : $core"
	echo
	echo "core is what is between \$PWD. and .pre_ltr.fasta"
	echo
	echo "Note: all files will be based on the working directory's name"
	echo
	exit
}


while [ $# -ne 0 ] ; do
	case $1 in
		-q|--q*|-m|--m*)
			shift; mapq=$1; shift ;;
		-i|--i*)
			shift; index=$1; shift ;;
		-c|--c*)
			shift; core=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#       Basically, this is TRUE AND DO ...
[ $# -gt 0 ] && usage

quality=`printf "Q%02d" $mapq`

date=`date "+%Y%m%d%H%M%S"`
log_file=`basename $0`.$quality.$index.$date.out

{

	find $basedir -maxdepth 2 \
		-name \*.$core.both_ltr.bowtie2.$index.$quality.insertion_points.overlappers \
		-exec cat {} \; > $index.$quality.insertion_points.overlappers
#	depth usage breaks on ec2
#	-depth will produce errors like the following on ec2. Use -maxdepth instead.
#	find: paths must precede expression: 2
#	Usage: find [-H] [-L] [-P] [-Olevel] [-D help|tree|search|stat|rates|opt|exec] [path...] [expression]
#		-depth 2 -exec cat {} \; > $index.$quality.insertion_points.overlappers

	find $basedir -maxdepth 2 \
		-name \*.$core.both_ltr.bowtie2.$index.$quality.rc_insertion_points.rc_overlappers \
		-exec cat {} \; > $index.$quality.rc_insertion_points.rc_overlappers
#	depth usage breaks on ec2
#	-depth will produce errors like the following on ec2. Use -maxdepth instead.
#	find: paths must precede expression: 2
#	Usage: find [-H] [-L] [-P] [-Olevel] [-D help|tree|search|stat|rates|opt|exec] [path...] [expression]
#		-depth 2 -exec cat {} \; > $index.$quality.rc_insertion_points.rc_overlappers

	awk '{print $2}' $index.$quality.insertion_points.overlappers \
		| sort -u > $index.$quality.insertion_points.overlappers.sort.uniq
	awk '{print $2}' $index.$quality.rc_insertion_points.rc_overlappers \
		| sort -u > $index.$quality.rc_insertion_points.rc_overlappers.sort.uniq


	awk '{print $0":F"}' $index.$quality.insertion_points.overlappers.sort.uniq \
		> overlapper_reference.$index.$quality
	awk '{print $0":R"}' $index.$quality.rc_insertion_points.rc_overlappers.sort.uniq \
		>> overlapper_reference.$index.$quality

	sort overlapper_reference.$index.$quality > overlapper_reference.$index.$quality.sorted

} 1>> $log_file 2>&1

#	not sure what was happening here, but works fine for me on my work laptop?
#	TYPO!  missed an underscore.  empty variable name.  "1>> 2>&1" is bad. 
#	./extract_insertion_points_and_overlappers.sh: line 116: $logfile: ambiguous redirect
#	the {} logfile doesn't seem to work with bash
