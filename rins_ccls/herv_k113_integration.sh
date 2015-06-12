#!/bin/sh -x
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


#	This really needs to be run in the data directory.
#
#	Eventually, may want to pass number of cpus or threads so
#	execs can use the same number.





#	If passed 1 fast[aq], check for chimeric reads.
#	If passed 2 fast[aq], also check for anchors with paired read run.

if [ $# -gt 2 -o $# -eq 0 ]; then
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <fasta/fastq file(s)>"
	echo
	echo "Expecting sorted/synchronised reads if paired given."
	echo
	exit
fi

{
	echo "Starting at ..."
	date

	#if [ $# -eq 2 ] ; then
	#	ln -s $1 raw.1.fastq
	#	ln -s $2 raw.2.fastq
	#fi
	#	else ...
	#	I expect the existance of raw.1.fastq and raw.2.fastq
	#

	#indexes=/Volumes/cube/working/indexes
	#	leading with the ": " stops execution
	#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
	#	would try to execute the result.  I just want the OR/EQUALS feature
	: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}


	#	used with all of the ec_fasta_split_and_blast calls
	#srun="srun --nice --share --exclude=n0000,n0001,n0002 --cpus-per-task=4"
	base=`basename $PWD`


	#	they MUST be exported, apparently, to be picked up by bowtie2
	export BOWTIE2_INDEXES


	#bowtie2
	#-x <bt2-idx> The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.


	[[ ${1:(-1)} == 'q' ]] && filetype='-q' || filetype='-f'

	if [ $# -eq 1 ] ; then
		files="-U $1"
	else
		files="-U $1,$2"
	fi

	base="$base.bowtie2.herv_k113_ltr_ends.__very_sensitive_local"

	bowtie2 --very-sensitive-local --threads 8 -x herv_k113_ltr_ends \
		$filetype $files -S $base.sam

	samtools view -b -S -F 4 -o $base.aligned.bam $base.sam
	rm $base.sam

	base="$base.aligned"

	samtools_extract_and_clip_chimeric_reads.sh $base.bam
	#	-> pre_ltr.fasta
	#	-> post_ltr.fasta

	bowtie2 -x hg19 --threads 8 -f $base.pre_ltr.fasta \
		-S $base.pre_ltr.bowtie2.hg19.sam
	#rm $base.pre_ltr.fasta

	samtools view -S -F 4 -b -o $base.pre_ltr.bowtie2.hg19.unsorted.bam \
		$base.pre_ltr.bowtie2.hg19.sam
	rm $base.pre_ltr.bowtie2.hg19.sam

	samtools sort $base.pre_ltr.bowtie2.hg19.unsorted.bam $base.pre_ltr.bowtie2.hg19
	rm $base.pre_ltr.bowtie2.hg19.unsorted.bam
	samtools index $base.pre_ltr.bowtie2.hg19.bam

	bowtie2 -x hg19 --threads 8 -f $base.post_ltr.fasta \
		-S $base.post_ltr.bowtie2.hg19.sam
	#rm $base.post_ltr.fasta

	samtools view -S -F 4 -b -o $base.post_ltr.bowtie2.hg19.unsorted.bam \
		$base.post_ltr.bowtie2.hg19.sam
	rm $base.post_ltr.bowtie2.hg19.sam

	samtools sort $base.post_ltr.bowtie2.hg19.unsorted.bam $base.post_ltr.bowtie2.hg19
	rm $base.post_ltr.bowtie2.hg19.unsorted.bam
	samtools index $base.post_ltr.bowtie2.hg19.bam

	#	find insertion points
	#	then find those with the signature overlap

	#	 4 = not aligned
	#	16 = reverse complement

	echo "Seeking insertion points and overlaps"

	samtools view -F 20 $base.pre_ltr.bowtie2.hg19.bam | awk '{print $3":"$4+length($10)}' | sort > $base.pre_ltr.bowtie2.hg19.insertion_points

	samtools view -F 20 $base.post_ltr.bowtie2.hg19.bam | awk '{print $3":"$4}' | sort > $base.post_ltr.bowtie2.hg19.insertion_points

	positions_within_10bp.sh $base.*.bowtie2.hg19.insertion_points | sort | uniq -c > $base.both_ltr.bowtie2.hg19.insertion_points.overlappers

	samtools view -F 4 -f 16 $base.pre_ltr.bowtie2.hg19.bam | awk '{print $3":"$4}' | sort > $base.pre_ltr.bowtie2.hg19.rc_insertion_points

	samtools view -F 4 -f 16 $base.post_ltr.bowtie2.hg19.bam | awk '{print $3":"$4+length($10)}' | sort > $base.post_ltr.bowtie2.hg19.rc_insertion_points

	positions_within_10bp.sh $base.*.bowtie2.hg19.rc_insertion_points | sort | uniq -c > $base.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers

	if [ $# -eq 2 ] ; then

		echo "Seeking paired-end anchors."

		base=`basename $PWD`
		base="$base.bowtie2.herv_k113"
		bowtie2 --threads 8 -x herv_k113 $filetype -1 $1 -2 $2 -S $base.sam

		#    f4F8 - Unmapped read whose mate did map
		samtools view -S -f 4 -F 8 -b -o $base.unaligned.bam $base.sam
		rm $base.sam

		samtools bam2fq $base.unaligned.bam > $base.unaligned.fastq
		rm $base.unaligned.bam

		bowtie2 --threads 8 -x hg19 $filetype $base.unaligned.fastq \
			-S $base.unaligned.bowtie2.hg19.sam
		rm $base.unaligned.fastq

		samtools view -S -b -F 4 -o $base.unaligned.bowtie2.hg19.aligned.bam \
			$base.unaligned.bowtie2.hg19.sam
		rm $base.unaligned.bowtie2.hg19.sam

	else
		echo "Only given 1 input file so not seeking paired-end anchors."
	fi

	base=`basename $PWD`
	samtools merge $base.herv_k113.hg19.aligned.unsorted.bam $base*hg19*bam
	samtools sort $base.herv_k113.hg19.aligned.unsorted.bam $base.herv_k113.hg19.aligned
	rm $base.herv_k113.hg19.aligned.unsorted.bam
	samtools index $base.herv_k113.hg19.aligned.bam

	echo
	echo "Finished at ..."
	date

} 1>>herv_k113_integration.out 2>&1
