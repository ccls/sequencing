#!/usr/bin/env bash

##!/bin/sh -x

index='hg19'
mapq=0

#	find pattern of base files NOT 
#		bam replaced with .pre_ltr.bowtie2.hg19.bam
#			            and .post_ltr.bowtie2.hg19.bam
findpattern='*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam'
basedir=`pwd`

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0`"
	echo
	echo "`basename $0` [--index STRING] [--mapq INTEGER] [--findpattern 'STRING']"
	echo
	echo "Example:"
	echo "  `basename $0`"
	echo
	echo "Defaults:"
	echo "  --index/-i ........ : $index"
	echo "  --mapq/-q  ........ : $mapq"
	echo "  --findpattern/-f .. : $findpattern"
	echo
	echo "Output directed to file, similarly named."
	echo
	exit 1
}

while [ $# -ne 0 ] ; do
	case $1 in
		-i|--i*)
			shift; index=$1; shift ;;
		-q|--mapq*)
			shift; mapq=$1; shift ;;
		-f|--f*)
			shift; findpattern=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

quality=`printf "Q%02d" $mapq`

date=`date "+%Y%m%d%H%M%S"`
logfile=`basename $0`.$index.$quality.$date.out

{
	echo "Starting at ..."
	date

	echo
	echo "index : $index"
	echo "mapq : $mapq"
	echo "quality : $quality"
	echo "findpattern : $findpattern"
	echo

	#TCGA-41-5651-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

	for bam in `find . -depth 2 -name $findpattern`; do
		cd $basedir

		echo $bam
		#	./TCGA-02-2483-01A/TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

		dirname $bam
		#	./TCGA-02-2483-01A
		cd `dirname $bam`

		pwd
		#	/Volumes/box/working/output/TCGA_Glioma_HERV52/TCGA-02-2483-01A

		base=`basename $bam`
		echo $base
		#	TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

		base=${base%.*}		#	drop the extension
		echo $base
		#	TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned


		#	TCGA-41-5651-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.bam

		samtools view -q $mapq -F 20 $base.pre_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4+length($10)}' \
			| sort > $base.pre_ltr.bowtie2.$index.$quality.insertion_points
		samtools view -q $mapq -F 20 $base.post_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4}' \
			| sort > $base.post_ltr.bowtie2.$index.$quality.insertion_points
		positions_within_10bp.sh $base.*.bowtie2.$index.$quality.insertion_points \
			| sort | uniq -c > $base.both_ltr.bowtie2.$index.$quality.insertion_points.overlappers

		samtools view -q $mapq -F 4 -f 16 $base.pre_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4}' \
			| sort > $base.pre_ltr.bowtie2.$index.$quality.rc_insertion_points
		samtools view -q $mapq -F 4 -f 16 $base.post_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4+length($10)}' \
			| sort > $base.post_ltr.bowtie2.$index.$quality.rc_insertion_points
		positions_within_10bp.sh $base.*.bowtie2.$index.$quality.rc_insertion_points \
			| sort | uniq -c > $base.both_ltr.bowtie2.$index.$quality.rc_insertion_points.rc_overlappers

	done

	echo
	echo "Finished at ..."
	date

} 1>> $logfile 2>&1
