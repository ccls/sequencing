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

	srun --nice --share --partition=bigmem \
		--job-name="Trinity_${name}" \
		--cpus-per-task=8 \
		--error=$base.trinity.errors.`date "+%Y%m%d%H%M%S"` \
		--output=$base.trinity.output.`date "+%Y%m%d%H%M%S"` \
		Trinity --CPU 8 \
			--bypass_java_version_check \
			--normalize_reads \
			--seqType fa \
			--run_as_paired \
			--min_contig_length 100 \
			--max_memory 40G \
			--single $1 \
			--output $base.trinity_output.nobackup &

	#		--full_cleanup \
	#	/my/home/ccls/data/working/TCGA_Glioma/trinity_wrapper.sh $1 &
	#mv $PWD/$base.trinity_output.nobackup.Trinity.fasta $PWD/$base.trinity.fasta

	#	Potential memory helpers
	#	--min_kmer_cov 2 \
	#	--normalize_reads

	shift
done