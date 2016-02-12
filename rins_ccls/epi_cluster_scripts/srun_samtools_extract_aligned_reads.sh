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
	name=${base##*/}	#	just in case given path

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	#if [[ ${ext,,} =~ sam ]] ; then
	#	flag='-S'
	#	samheader="-h $1"
	#else
	#	echo "Input file is not a sam file so not including header in output."
	#	flag=''
	#	samheader=''
	#fi
	[[ ${ext,,} =~ sam ]] && flag='-S' || flag=''

#	cmd="samtools view $flag -@ 8 -h -F 4 -b -o $base.aligned.bam $1"
#	echo $cmd

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="extract_aligned_${name}" \
		--cpus-per-task=8 \
		--error=$base.samtools_extract_aligned.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.samtools_extract_aligned.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		samtools_extract_aligned_reads.sh $1 &

	shift
done
