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

date=`date "+%Y%m%d%H%M%S"`

{
	echo "Starting at ..."
	date



#	Replace this for loop with a while loop reading the queue

	for line in `cat $1` ; do
		#	line ~ data/NA20505/sequence_read/ERR005686




		echo $line
		subject=${line#*/}
		subject=${subject%%/*}
		echo $subject
		sample=${line##*/}
		echo $sample

		mkdir $sample
		cd $sample

		#	fastq stuff seems to be under phase3/
		date
		aws s3 ls s3://1000genomes/phase3/${line}_1.filt.fastq.gz
		aws s3 ls s3://1000genomes/phase3/${line}_2.filt.fastq.gz
		aws s3 cp s3://1000genomes/phase3/${line}_1.filt.fastq.gz ./
		aws s3 cp s3://1000genomes/phase3/${line}_2.filt.fastq.gz ./
		date
		ls -trail
		gunzip *gz
		ls -trail
		date
		aws_fastq_to_herv_k113_overlappers.sh *fastq
		date
		ls -trail

		mkdir $sample
		mv *hg19.bam *insertion_points *overlappers *.out $sample/
		tar cfvz $sample.tar.gz $sample
		aws s3 cp $sample.tar.gz s3://sequers/1000genomes/$subject/

		cd ..
		/bin/rm -rf $sample

	done

	echo
	echo "Finished at ..."
	date
} 1>>$HOME/`basename $0`.$date.out 2>&1

aws s3 cp $HOME/`basename $0`.$date.out s3://sequers/1000genomes/

#	sudo shutdown -h now

