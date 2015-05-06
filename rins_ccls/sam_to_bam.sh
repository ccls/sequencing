#!/bin/sh

while [ $# -ne 0 ] ; do
	if [ -f "$1" ] ; then
		echo "Creating bam from (assumed) sam file $1"

		base=${1%.*}		#	drop the extension
		chmod -w $1

		echo "Converting to bam"
		#	Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
		samtools view -S -b -o $base.unsorted.bam $1
		chmod -w $base.unsorted.bam

		echo "Sorting"
		#	Usage: samtools sort [options...] [in.bam]
		#	Legacy usage: samtools sort [options...] <in.bam> <out.prefix>
		samtools sort $base.unsorted.bam $base
		chmod -w $base.bam

		echo "Indexing"
		#	Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
		samtools index $base.bam
		chmod -w $base.bam.bai
	else
		echo "$1 not a file? Skipping"
	fi
	shift
done
