#!/bin/bash

date=`date "+%Y%m%d%H%M%S"`
common="bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned"
dir="/Volumes/box/1000genomes/sequers/1000genomes/untarred/grouping"
cd $dir

prev_pop=""
group_size=3

for sample_w_pop in `cat select_samples_with_population` ; do

	echo $sample_w_pop
	sample=${sample_w_pop%%,*}
#	echo $sample
	pop=${sample_w_pop##*,}
#	echo $pop
	
	if [ -z $prev_pop ] ; then
		prev_pop=$pop
	fi
	if [ $pop != $prev_pop  ] ; then
		prev_pop=$pop
		group_list=( )
	fi

	#	add sample to group.  Array stuff needs curly braces.
	group_list=( ${group_list[*]} $sample )
#	echo ${group_list[*]}		#	entire array's contents
#	echo ${#group_list[*]}	#	array length

	if [[ ${#group_list[*]} -ge $group_size ]] ; then

#		echo "group found"
#		echo ${group_list[*]}		#	entire array's contents
		samples=${group_list[*]}	#	entire array's contents
		samples=${samples// /,}
#		echo $samples
		outname="$pop ${group_list[*]}"
		outname=${outname// /.}
		echo $outname

		#	make sure this name doesn't match the above find
		outdir="grouping.$date/$group_with_pop"
#		mkdir -p $outdir
	
		for ltr in pre_ltr post_ltr ; do
echo			ls -l {${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.bam
	
			if [ `ls -1 {${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.bam | wc -l` -gt 1 ] ; then
echo				samtools merge \
					$outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam \
					{${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.bam
			else
echo				cp {${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.bam \
					$outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam
			fi
echo			samtools sort \
				$outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam \
				$outdir/$outname.$common.$ltr.bowtie2.hg19	#	NOTE no .bam suffix
echo			\rm $outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam
echo			samtools index $outdir/$outname.$common.$ltr.bowtie2.hg19.bam
	
			for q in ALL Q10 Q20 ; do
				for ipt in insertion_points rc_insertion_points ; do
	
echo					ls -l {${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.$q.$ipt
echo					cat {${samples}}/{${samples}}.$common.$ltr.bowtie2.hg19.$q.$ipt \
						| sort > $outdir/$outname.$common.$ltr.bowtie2.hg19.$q.$ipt
	
				done	#	for ipt
			done	#	for q
		done	#	for ltr
	
		for q in ALL Q10 Q20 ; do
	
echo			ls -l $outdir/$outname.$common.*.bowtie2.hg19.$q.insertion_points
echo			positions_within_10bp.sh $outdir/$outname.$common.*.bowtie2.hg19.$q.insertion_points \
				| sort | uniq -c > $outdir/$outname.$common.both_ltr.bowtie2.hg19.$q.insertion_points.overlappers
	
			#	yes, expecting rc_overlappers
echo			ls -l $outdir/$outname.$common.*.bowtie2.hg19.$q.rc_insertion_points
echo			positions_within_10bp.sh $outdir/$outname.$common.*.bowtie2.hg19.$q.rc_insertion_points \
				| sort | uniq -c > $outdir/$outname.$common.both_ltr.bowtie2.hg19.$q.rc_insertion_points.rc_overlappers
	
		done	#	for q

		group_list=( )

	fi
done	#	for sample
