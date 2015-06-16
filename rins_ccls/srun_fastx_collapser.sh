#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}
	name=${base##*/}

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \
	srun --share --nice \
		--partition=bigmem \
		--exclude=n[0000-0009] \
		--job-name="fastx_collapser_${name}" \
		--output=$base.fastx_collapser.output.`date "+%Y%m%d%H%M%S"`.nobackup  \
		--error=$base.fastx_collapser.errors.`date "+%Y%m%d%H%M%S"`.nobackup  \
		fastx_collapser -i $1 -o $base.collapsed.fasta &


#	Not multithreaded so only blocks other jobs' access to node.
#		--cpus-per-task=8 \

	shift
done
