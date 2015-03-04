#!/bin/sh

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .fasta extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --share --nice --partition=bigmem \
		--begin=23:00 \
		--job-name="blastn_nt_${name}" \
		--cpus-per-task=8 \
		--output=$PWD/$base.blastn_nt.output.`date "+%Y%m%d%H%M%S"`  \
		--error=$PWD/$base.blastn_nt.errors.`date "+%Y%m%d%H%M%S"`  \
		blastn -num_threads 8 -num_alignments 20 -num_descriptions 20 \
			-evalue 0.05 -outfmt 0 -db nt \
			-query $1 \
			-out $base.blastn_nt.txt &

	shift
done
