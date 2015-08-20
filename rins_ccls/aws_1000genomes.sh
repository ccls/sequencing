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


if [ $# -ne 0 ]; then
	echo
	echo "Usage:"
	echo
	echo "`basename $0`"
	echo
	exit
fi

date=`date "+%Y%m%d%H%M%S"`

pid_file=$HOME/`basename $0`.$date.pid
echo $$ > $pid_file
log_file=$HOME/`basename $0`.$date.out

{
	echo "Starting at ..."
	date

	while [ -f $pid_file ]; do

		queue_url="https://us-west-1.queue.amazonaws.com/156714443422/1000genomes"

		message=`aws sqs receive-message --queue-url $queue_url`

		if [ -z "$message" ] ; then
			echo "Received message was blank.  I'm done."
			break	#	from the while loop
		fi

		line=`echo $message | python -c \
			'import sys, json; print json.load(sys.stdin)["Messages"][0]["Body"]'`
		echo $line

		handle=`echo $message | python -c \
			'import sys, json; print json.load(sys.stdin)["Messages"][0]["ReceiptHandle"]'`
		echo $handle

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

		#	apparently there is a limit on the number of inflight messages
		#	so delete them as soon as possible.
		aws sqs delete-message --queue-url $queue_url --receipt-handle $handle

	done	#	while [ -f $pid_file ]; do

	echo
	echo "Finished at ..."
	date
} 1>>$log_file 2>&1

aws s3 cp $log_file s3://sequers/1000genomes/

\rm $pid_file

#	sudo shutdown -h now

