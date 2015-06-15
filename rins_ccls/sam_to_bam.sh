#!/usr/bin/env bash

while [ $# -ne 0 ] ; do
	if [ -f "$1" ] ; then
		echo "Creating bam from (assumed) sam file $1"

		base=${1%.*}		#	drop the extension
		ext=${1##*.}		#	just the extension
#		chmod -w $1
		mv $1 $1.input
		chmod +w $base.bam
		chmod +w $base.bam.bai

		[[ ${ext,,} =~ sam ]] && flag='-S' || flag=''

		echo "Converting to bam"
		#	Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
		samtools view $flag -b -o $base.unsorted.bam $1.input
#		chmod -w $base.unsorted.bam
		if [ -s $base.unsorted.bam ] ; then rm $1.input ; fi

		echo "Sorting"
		#	Usage: samtools sort [options...] [in.bam]
		#	Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
		samtools sort -m 2G $base.unsorted.bam $base
#		chmod -w $base.bam
		if [ -s $base.bam ] ; then rm $base.unsorted.bam ; fi

		echo "Indexing"
		#	Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
		samtools index $base.bam

		chmod -w $base.bam
		chmod -w $base.bam.bai
	else
		echo "$1 not a file? Skipping"
	fi
	shift
done
