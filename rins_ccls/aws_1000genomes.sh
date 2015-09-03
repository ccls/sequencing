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

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--shutdown] [--threads INTEGER]"
	echo
	echo "--threads  : passed to bowtie2"
	echo "--shutdown : shutdown after complete"
	echo
	echo "Defaults:"
	echo "  threads ..... : 2"
	echo
	exit
}

threads=2
#shutdown='false'

date=`date "+%Y%m%d%H%M%S"`
pid_file=$HOME/`basename $0`.$date.pid
echo $$ > $pid_file
log_file=$HOME/`basename $0`.$date.`hostname`.$$.out
shutdown_file=$HOME/`basename $0`.shutdown

while [ $# -ne 0 ] ; do
	case $1 in
		-t|--t*)
			shift; threads=$1; shift;;
		-s|--s*)
			shift; touch $shutdown_file;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#       Basically, this is TRUE AND DO ...
[ $# -ne 0 ] && usage


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
		sample=${line#*/}
		sample=${sample%%/*}
		echo $sample
		submission=${line##*/}
		echo $submission

		#	aws_fastq_to_herv_k113_overlappers.sh is expecting 
		#	`basename $PWD` to be $submission
		mkdir -p working/$submission
		cd working/$submission

		#	fastq stuff seems to be under phase3/
		date
		aws s3 ls s3://1000genomes/phase3/${line}_1.filt.fastq.gz
		aws s3 ls s3://1000genomes/phase3/${line}_2.filt.fastq.gz
		aws s3 cp s3://1000genomes/phase3/${line}_1.filt.fastq.gz ./${sample}-${submission}_1.fastq.gz
		aws s3 cp s3://1000genomes/phase3/${line}_2.filt.fastq.gz ./${sample}-${submission}_2.fastq.gz
		date
		ls -trail
		gunzip *gz
		ls -trail
		date
		aws_fastq_to_herv_k113_overlappers.sh --threads $threads *fastq
		date
		ls -trail

		mkdir ${sample}-${submission}
		mv *hg19.bam *insertion_points *overlappers *.out ${sample}-${submission}/
		tar cfvz ${sample}-${submission}.tar.gz ${sample}-${submission}



#fdate=`date "+%Y%m%d%H%M%S"`
#		aws s3 ls s3://sequers/1000genomes/${sample}-${submission}.tar.gz

		aws s3 cp ${sample}-${submission}.tar.gz s3://sequers/1000genomes/





		cd ../..
		/bin/rm -rf working/$submission

#
#	I should really do some type of check to ensure it finished before deleting
#

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

#	sudo raises this error
#	sudo: sorry, you must have a tty to run sudo
#	adding -t, -tt, -ttt, -tttt to ssh doesn't change this but can get this instead....
#	Pseudo-terminal will not be allocated because stdin is not a terminal.
#	Only ...
# sudo visudo 
#   to comment out the following lines ...
# Defaults    requiretty
# Defaults   !visiblepw
#	works.  Need to create yet another AMI with this change.
#	http://unix.stackexchange.com/questions/49077
#if [ -f $shutdown_file -o $shutdown == 'true' ]; then
if [ -f $shutdown_file ]; then
	sudo shutdown -h now
	#	sudo halt
fi

