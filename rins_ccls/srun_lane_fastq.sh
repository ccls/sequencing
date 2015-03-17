#!/bin/sh

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .bam extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--exclude=n[0000-0009] \
		--job-name="lane_fastq_${name}" \
		--cpus-per-task=8 \
		--error=$base.lane_fastq.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.lane_fastq.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		awk '
			( /\/1$/ ){ f=sprintf("%s.1.fastq",substr(FILENAME,0,index(FILENAME,".fastq")-1)) }
			( /\/2$/ ){ f=sprintf("%s.2.fastq",substr(FILENAME,0,index(FILENAME,".fastq")-1)) }
			{ print >> f }
		' $1 &

	shift
done
