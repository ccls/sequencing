#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <bamfile(s) or samfile(s)>"
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
	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--partition=bigmem \
		--exclude=n[0000-0009] \
		--job-name="samtools_bam_to_paired_fastq_${name}" \
		--cpus-per-task=8 \
		--error=$base.samtools_bam_to_paired_fastq.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.samtools_bam_to_paired_fastq.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		samtools_bam_to_paired_fastq.sh $1 &

	shift
done
