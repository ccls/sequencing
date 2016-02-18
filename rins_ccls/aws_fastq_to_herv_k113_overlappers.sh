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


#	This really needs to be run in the data directory.
#
#	Eventually, may want to pass number of cpus or threads so
#	execs can use the same number.



#	If passed 1 fast[aq], check for chimeric reads.
#	If passed 2 fast[aq], also check for anchors with paired read run.

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--threads 4] <fastq file(s)>"
	echo
	echo "Defaults:"
	echo "  threads ..... : 2"
	echo
	echo "Expecting sorted/synchronised reads."
	echo
	echo "Note: all files will be based on the working directory's name"
	echo
	exit
}

threads=2

while [ $# -ne 0 ] ; do
	case $1 in
		-t|--t*)
			shift; threads=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done


#       Basically, this is TRUE AND DO ...
[ $# -gt 2 -o $# -eq 0 ] && usage



#	if no fastq files exist, base becomes "*fastq"
#base=`basename $1`
#base=${base%%_1.*}

base=`basename $PWD`

{
	echo "Starting at ..."
	date




#	wget bam file

#	sort bam? can i check if bam is sorted by name?  in header?
#
#	YES!
#
#	samtools view -H BAMFILE | grep @HD | 
#	
#	Tag Description @HD
#	The header line.  The first line if present.
#	VN * Format version.  Accepted format : /^[0-9]+\.[0-9]+$/ .
#	SO Sorting  order  of  alignments.  Valid  values : unknown (default), unsorted , queryname and coordinate .  For coordinate sort,  the major sort key is the RNAME field,  with order defined by the order of @SQ lines in the header.  The minor sort key is the POS field.  For alignments with equal RNAME and POS , order is arbitrary.  All alignments with ` * ' in RNAME field follow alignments with some other value but otherwise are in arbitrary order.
#	GO Grouping of alignments, indicating that similar alignment records are grouped together but the file is not necessarily sorted overall.  Valid values : none (default), query (alignments are grouped by QNAME), and reference (alignments are grouped by RNAME / POS).
#	@HD	VN:1.0	SO:queryname
#	bam_to_paired_fastq


	#indexes=/Volumes/cube/working/indexes
	#	leading with the ": " stops execution
	#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
	#	would try to execute the result.  I just want the OR/EQUALS feature
	#: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}
	: ${BOWTIE2_INDEXES:="$HOME/BOWTIE2_INDEXES"}


	#	used with all of the ec_fasta_split_and_blast calls
	#srun="srun --nice --share --exclude=n0000,n0001,n0002 --cpus-per-task=4"


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

	

#	bam file - extension + _R?.fastq
#
#	echo $1
#	core=${1%.*}		#	drop the extension
#	files="-U ${core}_R1.fastq,${core}_R2.fastq"
#	filetype='-q'
#



	base="$base.bowtie2.herv_k113_ltr_ends.__very_sensitive_local"
	bowtie2 --very-sensitive-local --threads $threads -x herv_k113_ltr_ends \
		$filetype $files -S $base.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	#	For the HUGE files, delete the fastq files once the sam file is done.
	rm $1
	rm $2

	samtools view -b -S -F 4 -o $base.aligned.bam $base.sam
	rm $base.sam

	base="$base.aligned"

	samtools_extract_and_clip_chimeric_reads.sh $base.bam
	#	-> pre_ltr.fasta
	#	-> post_ltr.fasta

	bowtie2 -x hg19 --threads $threads -f $base.pre_ltr.fasta \
		-S $base.pre_ltr.bowtie2.hg19.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	bowtie2 -x hg19 --threads $threads -f $base.post_ltr.fasta \
		-S $base.post_ltr.bowtie2.hg19.sam
	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	#	find insertion points
	#	then find those with the signature overlap

	#	 f = ALL/YES
	#	 F = NONE/NOT	(results in double negatives)
	#	 4 = not aligned
	#	 8 = mate not aligned
	#	16 = reverse complement

	echo "Seeking insertion points and overlaps"

	samtools view -F 20 $base.pre_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.pre_ltr.bowtie2.hg19.Q00.insertion_points
	samtools view -F 20 $base.post_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4}' \
		| sort > $base.post_ltr.bowtie2.hg19.Q00.insertion_points
	positions_within_10bp.sh $base.*.bowtie2.hg19.Q00.insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.hg19.Q00.insertion_points.overlappers

	samtools view -F 4 -f 16 $base.pre_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4}' \
		| sort > $base.pre_ltr.bowtie2.hg19.Q00.rc_insertion_points
	samtools view -F 4 -f 16 $base.post_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.post_ltr.bowtie2.hg19.Q00.rc_insertion_points
	positions_within_10bp.sh $base.*.bowtie2.hg19.Q00.rc_insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.hg19.Q00.rc_insertion_points.rc_overlappers


	samtools view -q 20 -F 20 $base.pre_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.pre_ltr.bowtie2.hg19.Q20.insertion_points
	samtools view -q 20 -F 20 $base.post_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4}' \
		| sort > $base.post_ltr.bowtie2.hg19.Q20.insertion_points
	positions_within_10bp.sh $base.*.bowtie2.hg19.Q20.insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.hg19.Q20.insertion_points.overlappers

	samtools view -q 20 -F 4 -f 16 $base.pre_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4}' \
		| sort > $base.pre_ltr.bowtie2.hg19.Q20.rc_insertion_points
	samtools view -q 20 -F 4 -f 16 $base.post_ltr.bowtie2.hg19.sam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.post_ltr.bowtie2.hg19.Q20.rc_insertion_points
	positions_within_10bp.sh $base.*.bowtie2.hg19.Q20.rc_insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.hg19.Q20.rc_insertion_points.rc_overlappers


	samtools view -S -b -o $base.pre_ltr.bowtie2.hg19.bam  $base.pre_ltr.bowtie2.hg19.sam
	rm $base.pre_ltr.bowtie2.hg19.sam

	samtools view -S -b -o $base.post_ltr.bowtie2.hg19.bam $base.post_ltr.bowtie2.hg19.sam
	rm $base.post_ltr.bowtie2.hg19.sam

	echo
	echo "Finished at ..."
	date

} 1>>$base.`basename $0`.out 2>&1

