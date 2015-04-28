#!/bin/sh

while [ $# -ne 0 ] ; do
	if [ -f "$1" ] ; then
		echo "Creating bam from (assumed) sam file $1"

		base=${1%.*}		#	drop the extension

		echo "Converting to bam"
		samtools view -S -b -o $base.unsorted.bam $1

		echo "Sorting"
		samtools sort $base.unsorted.bam $base.bam

		echo "Indexing"
		samtools index $base.bam
	else
		echo "$1 not a file? Skipping"
	fi
	shift
done
