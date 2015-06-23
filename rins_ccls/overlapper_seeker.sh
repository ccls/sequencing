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
	: ${BOWTIE2_INDEXES:="/Volumes/box/indexes"}


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

	for fastafile in $BOWTIE2_INDEXES/overlappers/*.fasta ; do
		echo $fastafile
		fastafile=${fastafile##*/}
		echo $fastafile
		fastabase=${fastafile%%.*}
		echo $fastabase

#		continue
	
		base="$base.bowtie2.$fastabase.__very_sensitive_local"
		bowtie2 --very-sensitive-local --threads 8 -x overlappers/$fastabase \
			$filetype $files -S $base.sam
		samtools view -b -S -F 4 -o $base.aligned.unsorted.bam $base.sam
		rm $base.sam
		samtools sort $base.aligned.unsorted.bam $base.aligned
		rm $base.aligned.unsorted.bam
		samtools index $base.aligned.bam
	
		base="$base.aligned"
	
		#	--posix NEEDS to be AFTER any -v settings!
		samtools view $flag -h -F 4 $base.bam | awk -v base=$base --posix '
			BEGIN {
				pre_out=sprintf("%s.pre_ltr.fasta",base)
				post_out=sprintf("%s.post_ltr.fasta",base)
			}
			( ( NR % 10000 ) == 0 ){ print "Read "NR" records" }
	
			( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4) }
	
			( ( $6 ~ /^[0-9]{2,}S[0-9IDM]*$/ ) && ( $4 <= 5 ) ){
				split($6,a,"S")
				clip=a[1]-$4+1
				print ">"$1"_pre_ltr" >> pre_out
				print substr($10,1,clip) >> pre_out
			}
			( ( $6 ~ /^[0-9IDM]*[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - (length($10)*0.9) ) ) ){
				clip=ref[$3]-$4+2
				print ">"$1"_post_ltr" >> post_out
				print substr($10,clip) >> post_out
			}'
		#	-> pre_ltr.fasta
		#	-> post_ltr.fasta
	
		if [ -f $base.pre_ltr.fasta ] ; then
			bowtie2 -x hg19 --threads 8 -f $base.pre_ltr.fasta \
				-S $base.pre_ltr.bowtie2.hg19.sam
			#rm $base.pre_ltr.fasta
			samtools view -S -F 4 -b -o $base.pre_ltr.bowtie2.hg19.unsorted.bam \
				$base.pre_ltr.bowtie2.hg19.sam
			rm $base.pre_ltr.bowtie2.hg19.sam
			samtools sort $base.pre_ltr.bowtie2.hg19.unsorted.bam $base.pre_ltr.bowtie2.hg19
			rm $base.pre_ltr.bowtie2.hg19.unsorted.bam
			samtools index $base.pre_ltr.bowtie2.hg19.bam
			samtools view -F 20 $base.pre_ltr.bowtie2.hg19.bam | \
				awk '{print $3":"$4+length($10)}' | \
				sort > $base.pre_ltr.bowtie2.hg19.insertion_points
			samtools view -F 4 -f 16 $base.pre_ltr.bowtie2.hg19.bam | \
				awk '{print $3":"$4}' | sort > $base.pre_ltr.bowtie2.hg19.rc_insertion_points
		fi
		
		if [ -f $base.post_ltr.fasta ] ; then
			bowtie2 -x hg19 --threads 8 -f $base.post_ltr.fasta \
				-S $base.post_ltr.bowtie2.hg19.sam
			#rm $base.post_ltr.fasta
			samtools view -S -F 4 -b -o $base.post_ltr.bowtie2.hg19.unsorted.bam \
				$base.post_ltr.bowtie2.hg19.sam
			rm $base.post_ltr.bowtie2.hg19.sam
			samtools sort $base.post_ltr.bowtie2.hg19.unsorted.bam $base.post_ltr.bowtie2.hg19
			rm $base.post_ltr.bowtie2.hg19.unsorted.bam
			samtools index $base.post_ltr.bowtie2.hg19.bam
			samtools view -F 20 $base.post_ltr.bowtie2.hg19.bam | \
				awk '{print $3":"$4}' | \
				sort > $base.post_ltr.bowtie2.hg19.insertion_points
			samtools view -F 4 -f 16 $base.post_ltr.bowtie2.hg19.bam | \
				awk '{print $3":"$4+length($10)}' | \
				sort > $base.post_ltr.bowtie2.hg19.rc_insertion_points
		fi
	
		#	find insertion points
		#	then find those with the signature overlap
	
		#	 f = ALL/YES
		#	 F = NONE/NOT	(results in double negatives)
		#	 4 = not aligned
		#	 8 = mate not aligned
		#	16 = reverse complement
	
		if [ -f $base.pre_ltr.bowtie2.hg19.insertion_points -a \
			-f $base.post_ltr.bowtie2.hg19.insertion_points ] ; then
			positions_within_10bp.sh $base.*.bowtie2.hg19.insertion_points | \
				sort | uniq -c > $base.both_ltr.bowtie2.hg19.insertion_points.overlappers
		fi
	
		if [ -f $base.pre_ltr.bowtie2.hg19.rc_insertion_points -a \
			-f $base.post_ltr.bowtie2.hg19.rc_insertion_points ] ; then
			positions_within_10bp.sh $base.*.bowtie2.hg19.rc_insertion_points | \
				sort | uniq -c > $base.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
		fi

#		#		longbase="$base.bowtie2.$fastabase.__very_sensitive_local.aligned"
#		longbase=$base
#		base=`basename $PWD`
#
#		if [ -f $longbase.pre_ltr.bowtie2.hg19.bam -a -f $longbase.post_ltr.bowtie2.hg19.bam ] ; then
#			samtools merge $base.$fastabase.hg19.aligned.unsorted.bam $longbase.*hg19*bam
#			samtools sort $base.$fastabase.hg19.aligned.unsorted.bam $base.$fastabase.hg19.aligned
#			rm $base.$fastabase.hg19.aligned.unsorted.bam
#			samtools index $base.$fastabase.hg19.aligned.bam
#		fi


#		#	Do this AFTER the merging, so it doesn't get merged (or be more specific)
#		samtools view -S -b -o $base.bowtie2.hg19.unsorted.bam $base.bowtie2.hg19.sam
#		rm $base.bowtie2.hg19.sam
#		samtools sort $base.bowtie2.hg19.unsorted.bam $base.bowtie2.hg19
#		rm $base.bowtie2.hg19.unsorted.bam
#		samtools index $base.bowtie2.hg19.bam

	done	#	for f in $BOWTIE2_INDEXES/overlappers/*.fasta ;

	echo
	echo "Finished at ..."
	date

} 1>>overlapper_seeker.out 2>&1
