#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <samfile(s) or bamfile(s)>"
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
#	name=${base##*/}	#	just in case given path

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	[[ ${ext,,} =~ sam ]] && flag='-S' || flag=''
	#echo $ext
	#echo $flag

	#    F4 = NOT unmapped = mapped
	samtools view $flag -@ 8 -b -F 4 -o $base.aligned.unsorted.bam $1
	samtools sort $base.aligned.unsorted.bam $base.aligned
	chmod -w $base.aligned.bam
	samtools index $base.aligned.bam
	chmod -w $base.aligned.bam.bai
	rm $base.aligned.unsorted.bam

	shift
done
