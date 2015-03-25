#!/bin/sh

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}
	name=${base#*/}

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \
	srun --share --nice \
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="fastq_to_fasta_${name}" \
		--output=$PWD/$base.fastq_to_fasta.output.`date "+%Y%m%d%H%M%S"`.nobackup  \
		--error=$PWD/$base.fastq_to_fasta.errors.`date "+%Y%m%d%H%M%S"`.nobackup  \
		fastq_to_fasta -Q33 -n -i $PWD/$1 -o $PWD/$base.fasta &

#	Not multithreaded so only blocks other jobs' access to node.
#		--cpus-per-task=8 \

	shift
done
