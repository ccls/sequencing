#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <samfile(s) or bamfile(s)>"
	echo
	echo "Example:"
	echo "  `basename $0` /my/path/*bam"
	echo
	echo "Assuming ..."
	echo "      ... input was aligned with bowtie2 and a form of --local"
	echo "          to a database with just 2 reads longer than the input reads."
	echo "          1 db entry is the 'beginning' of the start of the ltr."
	echo "          1 db entry is the 'ending' of the end of the ltr."
	echo
	echo "Output is currently just 2 fasta files containing the pre- and post-ltr chimeric reads."
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
	name=${base#*/}	#	just in case given path, drop the path

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	#if [[ ${ext,,} =~ sam ]] ; then
	#	flag='-S'
	#	samheader="-h $1"
	#else
	#	echo "Input file is not a sam file so not including header in output."
	#	flag=''
	#	samheader=''
	#fi
	[[ ${ext,,} =~ sam ]] && flag='-S' || flag=''
	#echo $ext
	#echo $flag


	#	The following numbers are all based on 100bp reads. (101bp actually)
	#	Every read.  Hmm.  Make it computed value as reads not always all the same length.

	#	
	#	Find alignments that align past the appropriate end of the ends of the ltr.
	#

	#    f4 = unmapped
	#    F4 = NOT unmapped = mapped
	#    F8 = mate NOT unmapped = mate mapped

	#    select LEFT soft clips from sam into fasta  ( ##S##M )
	#	Fastq?
#	samtools view $flag -h -F 4 $1 | awk '
#	( ( $6 ~ /^[0-9]{2}S[0-9]{1,2}M$/ ) && ( $4 <= 5 ) && ( $3 ~ /beginning/ ) ){
#		split($6,s,"S")
#		clip=s[1]-$4+1
#		print ">"$1”_pre_ltr”
#		print substr($10,1,clip)
#	}' > $base.pre_ltr.fasta
#
#
#	#	Could do both of these in 1 awk call .... just sayin'
#	#	BEGIN{ 
#	#	pre="$base.pre_ltr.fasta"
#	#	post="$base.post_ltr.fasta"
#
#	#	could also read reference lengths from sam header
#	#	@SQ	SN:HERV52I	LN:3556
#	#	@SQ	SN:HERV57I	LN:5343
#	#	@SQ	SN:HERV9	LN:8399
#	#	@SQ	SN:HERVE	LN:7813
#	#	@SQ	SN:HERVE_a	LN:7847
#	#	samtools view -H FILE | awk '( /^@SQ/ ){ sn[substr($2,4)] = substr($3,4) }END{ for( key in sn ){ print key,sn[key]}}'
#
##	ref[$3] = reference length
#
#
#	#    select RIGHT soft clips from sam into fasta  ( ##M##S )
#	samtools view $flag -h -F 4 $1 | awk '
#	function sum(a){
#		s=0;
#		for(i=0;i<length(a);i++){s+=a[i]};
#		return s;
#	}
#	( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4) }
#	( ( $6 ~ /^[0-9]{1,2}M[0-9]{2}S$/ ) && ( $4 >= 55 ) && ( $3 ~ /ending/ ) ){
#		split($6,a,"[[:alpha:]]")
#		clip=sum(a)-$4+1
#		print ">"$1”_post_ltr”
#		print substr($10,clip)
#	}' > $base.post_ltr.fasta

#	function sum(a){
#		s=0;
#		for(i=0;i<length(a);i++){s+=a[i]};
#		return s;
#	}

	samtools view $flag -h -F 4 $1 | awk '
	( /^@SQ/ ){ ref[substr($2,4)] = substr($3,4) }

	( ( $6 ~ /^[0-9]{2,}S[0-9]{1,}M$/ ) && ( $4 <= 5 ) && ( $3 ~ /beginning/ ) ){
		split($6,a,"S")
		clip=a[1]-$4+1
		print ">"$1"_pre_ltr" >> "'$base.pre_ltr.fasta'"
		print substr($10,1,clip) >> "'$base.pre_ltr.fasta'"
	}

	( ( $6 ~ /^[0-9]{1,}M[0-9]{2,}S$/ ) && ( $4 >= ( ref[$3] - (length($10)*0.8) ) && ( $3 ~ /ending/ ) )){
		clip=ref[$3]-$4+2
		print ">"$1"_post_ltr" >> "'$base.post_ltr.fasta'"
		print substr($10,clip) >> "'$base.post_ltr.fasta'"
	}'

#		split($6,a,"[[:alpha:]]")
#		#clip=sum(a)-$4+1




	#	bowtie2 -x hg19 -U $base.pre_ltr.fasta  -S $base.pre_ltr.bowtie2.hg19.sam
	#	bowtie2 -x hg19 -U $base.post_ltr.fasta -S $base.post_ltr.bowtie2.hg19.sam

#    $4+$6 correctly does integer math and ignores the trailing M from $6
#samtools view -q 2 TCGA-FG-7636-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.beginning.soft_clipped.bowtie2.hg19..bam | awk '{print $3,$4,$6,$4+$6}'

#samtools view -q 2 TCGA-FG-7636-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.ending.soft_clipped.bowtie2.hg19..bam | awk '{print $3,$4,$6,$4+$6}' > TCGA-FG-7636-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.ending.soft_clipped.bowtie2.hg19..insertion_points


	shift
done
