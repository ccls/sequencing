#!/bin/bash -x
##!/usr/bin/env bash -x
#
#	Can't pass options when using env
#
#	with the -x the commands are sent to STDERR before execution
#
#	bowtie output goes to stderr for some reason
#	probably because the SAM file usually goes to stdout
#	so, wrap everything in curly braces and direct both
#	to files.
#
#	Explicit redirection within the block will override this
#


if [ $# -ne 1 ]; then
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <list file>"
	echo
	echo "list file lines like ..."
	echo "data/NA20505/sequence_read/ERR005686"
	echo
	exit
fi


{
	echo "Starting at ..."
	date

	for line in `cat $1` ; do
		#	line ~ data/NA20505/sequence_read/ERR005686
		echo $line
		sample=${line##*/}
		echo $sample

		mkdir $sample
		cd $sample

	#	#	fastq stuff seems to be under phase3/
	#	aws s3 cp s3://1000genomes/phase3/${line}_1.filt.fastq.gz ./
	#	aws s3 cp s3://1000genomes/phase3/${line}_2.filt.fastq.gz ./
	#	gunzip *gz
	#	aws_fastq_to_herv_k113_overlappers.sh *fastq
	#	mkdir $sample
	#	mv *hg19.bam *insertion_points *overlappers *.out $sample/
	#	tar cfvz $sample.tar.gz $sample
	#	aws s3 cp $sample.tar.gz s3://sequers/1000genomes/

		cd ..
		/bin/rm -rf $sample

	done

	echo
	echo "Finished at ..."
	date
} #	1>>`basename $0`.out 2>&1
