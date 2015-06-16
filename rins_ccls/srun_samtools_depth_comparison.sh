#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` bamfile1 bamfile2"
	echo
	echo "Example:"
	echo "  `basename $0` bamfile1 bamfile2"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ge 2 ] ; do
	echo $1
	echo $2

	base=${1%.*}		#	drop the extension
#	ext=${1##*.}		#	grab the extension
	name=${base##*/}	#	just in case given path

	cmd="samtools_depth_comparison.sh $1 $2"
	echo $cmd

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="samtools_depth_comparison_${name}" \
		--cpus-per-task=8 \
		--error=$base.samtools_depth_comparison.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.samtools_depth_comparison.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		$cmd &

	shift
	shift
done
