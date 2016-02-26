#!/usr/bin/env bash
set -x


index='hg19'
mapq=0
core='bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned'
basedir=`pwd`


function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0`"
	echo
	echo "`basename $0` [--index STRING] [--mapq INTEGER] [--core 'STRING'] sampledir(s)"
	echo
	echo "Example:"
	echo "  `basename $0` NA211??"
	echo "  `basename $0` -i hg38 -q 30 NA2114?"
	echo
	echo "Defaults:"
	echo "  --index/-i .. : $index"
	echo "  --mapq/-q  .. : $mapq"
	echo "  --core/-c ... : $core"
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
		-c|--c*)
			shift; core=$1; shift ;;
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
	echo "core : $core"
	echo

	#TCGA-41-5651-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

	while [ $# -ne 0 ] ; do
		cd $basedir

		#echo $bam
		#	./TCGA-02-2483-01A/TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

		#dirname $bam
		#	./TCGA-02-2483-01A
		#cd `dirname $bam`
		cd $1

		pwd
		#	/Volumes/box/working/output/TCGA_Glioma_HERV52/TCGA-02-2483-01A

		#base=`basename $bam`
		#echo $base
		#	TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

		echo $core
		#	*bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned

		#base=${base%.*}		#	drop the extension
		#echo $base
		#	TCGA-02-2483-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned

		sample=`basename $PWD`


		#	TCGA-41-5651-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.bam

		samtools view -q $mapq -F 20 $sample.$core.pre_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4+length($10)}' \
			| sort > $sample.$core.pre_ltr.bowtie2.$index.$quality.insertion_points
		samtools view -q $mapq -F 20 $sample.$core.post_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4}' \
			| sort > $sample.$core.post_ltr.bowtie2.$index.$quality.insertion_points
		positions_within_10bp.sh $sample.$core.*.bowtie2.$index.$quality.insertion_points \
			| sort | uniq -c > $sample.$core.both_ltr.bowtie2.$index.$quality.insertion_points.overlappers

		samtools view -q $mapq -F 4 -f 16 $sample.$core.pre_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4}' \
			| sort > $sample.$core.pre_ltr.bowtie2.$index.$quality.rc_insertion_points
		samtools view -q $mapq -F 4 -f 16 $sample.$core.post_ltr.bowtie2.$index.bam \
			| awk '{print $3":"$4+length($10)}' \
			| sort > $sample.$core.post_ltr.bowtie2.$index.$quality.rc_insertion_points
		positions_within_10bp.sh $sample.$core.*.bowtie2.$index.$quality.rc_insertion_points \
			| sort | uniq -c > $sample.$core.both_ltr.bowtie2.$index.$quality.rc_insertion_points.rc_overlappers

		shift
	done

	echo
	echo "Finished at ..."
	date

} 1>> $logfile 2>&1
