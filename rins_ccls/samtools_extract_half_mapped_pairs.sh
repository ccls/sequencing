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
	name=${base##*/}	#	just in case given path

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	if [[ ${ext,,} =~ sam ]] ; then
		flag='-S'
		samheader="-h $1"
	else
		echo "Input file is not a sam file so not including header in output."
		flag=''
		samheader=''
	fi
	#[[ ${ext,,} =~ sam ]] && flag='-S' || flag=''
	#echo $ext
	#echo $flag

	#    F4 = NOT unmapped = mapped
	samtools view $flag -@ 8 -b -F 4 -o $name.mapped.bam $1

	#    f4 = unmapped
	#    F8 = mate NOT unmapped
	samtools view $flag -@ 8 -b -f 4 -F 8 -o $name.mappedmate.bam $1

	#[[ ${ext,,} =~ sam ]] && samheader="-h $1" || samheader=''

	#	older versions don't seem to actually sort
	samtools merge -n -@ 8 $samheader \
		$name.MERGED.bam \
		$name.mapped.bam \
		$name.mappedmate.bam

	samtools sort -n $name.MERGED.bam $name.MERGEDANDSORTED

#	samtools bam2fq $name.MERGEDANDSORTED.bam > $name.MERGEDANDSORTED.fastq

	shift
done
