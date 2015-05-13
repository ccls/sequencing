#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` bamfile(s)"
	echo
	echo "Example:"
	echo "  `basename $0` /my/path/*bam"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .bam extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--partition=bigmem \
		--exclude=n[0000-0009] \
		--job-name="samtools_sort_${name}" \
		--cpus-per-task=8 \
		--error=$base.samtools_sort.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.samtools_sort.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		samtools sort -n -m 3G -@ 8 $1 $base.sortedbyname &

	shift
done
