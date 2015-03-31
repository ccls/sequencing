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
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="lane_fastq_${name}" \
		--cpus-per-task=8 \
		--error=$base.lane_fastq.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.lane_fastq.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		awk '{
			if( ( (NR-1) % 4 ) == 0 ){
				f=sprintf("%s.%i.fastq", 
					substr(FILENAME,0,index(FILENAME,".fastq")-1), 
					substr($0,length,1))
			}
			print >> f 
		}' $1 &

	shift
done

#	sadly, fastq allows @ and / as well as 1 and 2 to be quality
#	scores so the naming scheme of starting with @ and ending with /1 or /2
#	as a defining characteristic is not enough.  Fastq is also
#	semi-defined as being 4 lines, mostly because the reads are
#	usually short.  Contigs will always be fasta and not fastq.
#	So the following won't work.
#
#awk '
#	( /\/1$/ ){ f=sprintf("%s.1.fastq",substr(FILENAME,0,index(FILENAME,".fastq")-1)) }
#	( /\/2$/ ){ f=sprintf("%s.2.fastq",substr(FILENAME,0,index(FILENAME,".fastq")-1)) }
#	{ print >> f }
#		' $1 &
#
#	So basically, just take the last character of the first line of 
#	every 4 line block to determine which lane it is.
