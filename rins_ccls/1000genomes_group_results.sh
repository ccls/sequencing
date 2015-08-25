#!/bin/bash

date=`date "+%Y%m%d%H%M%S"`

for sample in `find . -type d -name \*-\* -depth 1 -exec basename {} \; | awk -F- '{print $1}' | uniq`; do

#for sample in `find . -type d -name HG000\*-\* -exec basename {} \; | awk -F- '{print $1}' | uniq`; do
#for sample in HG00096 HG00097 HG00099 HG00100 HG00103 HG00104 HG00106 ; do
#for sample in HG00096 HG00097 HG00099 ; do

	echo $sample

	common="bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned"

	#	don't match the above find
	outdir="grouping.$date/$sample"
	mkdir -p $outdir

	for ltr in pre_ltr post_ltr ; do
		ls -l ${sample}-*/${sample}*.$common.$ltr.bowtie2.hg19.bam

		if [ `ls -1 ${sample}-*/${sample}*.$common.$ltr.bowtie2.hg19.bam | wc -l` -gt 1 ] ; then
			samtools merge \
				$outdir/$sample.$common.unsorted.$ltr.bowtie2.hg19.bam \
				${sample}-*/${sample}*.$common.$ltr.bowtie2.hg19.bam
		else
			cp ${sample}-*/${sample}*.$common.$ltr.bowtie2.hg19.bam \
				$outdir/$sample.$common.unsorted.$ltr.bowtie2.hg19.bam
		fi
		samtools sort \
			$outdir/$sample.$common.unsorted.$ltr.bowtie2.hg19.bam \
			$outdir/$sample.$common.$ltr.bowtie2.hg19	#	NOTE no .bam suffix
		\rm $outdir/$sample.$common.unsorted.$ltr.bowtie2.hg19.bam
		samtools index $outdir/$sample.$common.$ltr.bowtie2.hg19.bam

		for q in ALL Q20 ; do
			for ipt in insertion_points rc_insertion_points ; do

				cat ${sample}-*/${sample}-*.$common.$ltr.bowtie2.hg19.$q.$ipt \
					| sort > $outdir/$sample.$common.$ltr.bowtie2.hg19.$q.$ipt

			done	#	for ipt
		done	#	for q
	done	#	for ltr

	for q in ALL Q20 ; do

		positions_within_10bp.sh $outdir/$sample.$common.*.bowtie2.hg19.$q.insertion_points \
			| sort | uniq -c > $outdir/$sample.$common.both_ltr.bowtie2.hg19.$q.insertion_points.overlappers

		#	yes, expecting rc_overlappers
		positions_within_10bp.sh $outdir/$sample.$common.*.bowtie2.hg19.$q.rc_insertion_points \
			| sort | uniq -c > $outdir/$sample.$common.both_ltr.bowtie2.hg19.$q.rc_insertion_points.rc_overlappers

	done	#	for q
done	#	for sample
