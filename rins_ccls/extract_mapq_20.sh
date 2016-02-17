#!/bin/sh -x

index=hg19
basedir=`pwd`


{
echo "Starting at ..."
date


#TCGA-41-5651-01A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam

for bam in `find . -depth 2 -name '*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.bam'`; do
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

	samtools view -q 20 -F 20 $base.pre_ltr.bowtie2.$index.bam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.pre_ltr.bowtie2.$index.Q20.insertion_points
	samtools view -q 20 -F 20 $base.post_ltr.bowtie2.$index.bam \
		| awk '{print $3":"$4}' \
		| sort > $base.post_ltr.bowtie2.$index.Q20.insertion_points
	positions_within_10bp.sh $base.*.bowtie2.$index.Q20.insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.$index.Q20.insertion_points.overlappers

	samtools view -q 20 -F 4 -f 16 $base.pre_ltr.bowtie2.$index.bam \
		| awk '{print $3":"$4}' \
		| sort > $base.pre_ltr.bowtie2.$index.Q20.rc_insertion_points
	samtools view -q 20 -F 4 -f 16 $base.post_ltr.bowtie2.$index.bam \
		| awk '{print $3":"$4+length($10)}' \
		| sort > $base.post_ltr.bowtie2.$index.Q20.rc_insertion_points
	positions_within_10bp.sh $base.*.bowtie2.$index.Q20.rc_insertion_points \
		| sort | uniq -c > $base.both_ltr.bowtie2.$index.Q20.rc_insertion_points.rc_overlappers

done

echo
echo "Finished at ..."
date
} 1>> `basename $0`.out 2>&1
