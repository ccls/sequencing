#!/bin/bash

date=`date "+%Y%m%d%H%M%S"`
common="bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned"
dir="/Volumes/box/1000genomes/sequers/1000genomes/untarred/grouping"
cd $dir

prev_pop=""
group_size=3

for sample_w_pop in `cat /Volumes/box/1000genomes/select_samples_with_population` ; do

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
		outname=${outname// /-}
		echo $outname

		#	make sure this name doesn't match the above find
		outdir="grouping.$date/$outname"
		mkdir -p $outdir
	
		for ltr in pre_ltr post_ltr ; do
			echo ls -l {${samples}}/*.$common.$ltr.bowtie2.hg19.bam

			#	Can't easily do brace expansion inside script.  So ...
			bams=${group_list[*]/%//*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.$ltr.bowtie2.hg19.bam}
			ls -l $bams
	
			if [ `ls -1 $bams | wc -l` -gt 1 ] ; then
				samtools merge \
					$outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam $bams
			else
				cp $bams $outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam
			fi

			samtools sort \
				$outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam \
				$outdir/$outname.$common.$ltr.bowtie2.hg19	#	NOTE no .bam suffix
			\rm $outdir/$outname.$common.unsorted.$ltr.bowtie2.hg19.bam
			samtools index $outdir/$outname.$common.$ltr.bowtie2.hg19.bam
	
			for q in ALL Q10 Q20 ; do
				for ipt in insertion_points rc_insertion_points ; do
	
					ipts=${group_list[*]/%//*.$common.$ltr.bowtie2.hg19.$q.$ipt}
					ls -l $ipts
					cat $ipts | sort > $outdir/$outname.$common.$ltr.bowtie2.hg19.$q.$ipt
	
				done	#	for ipt
			done	#	for q
		done	#	for ltr

		for q in ALL Q10 Q20 ; do
	
			ls -l $outdir/$outname.$common.*.bowtie2.hg19.$q.insertion_points
			positions_within_10bp.sh $outdir/$outname.$common.*.bowtie2.hg19.$q.insertion_points \
				| sort | uniq -c > $outdir/$outname.$common.both_ltr.bowtie2.hg19.$q.insertion_points.overlappers

			#	yes, expecting rc_overlappers
			ls -l $outdir/$outname.$common.*.bowtie2.hg19.$q.rc_insertion_points
			positions_within_10bp.sh $outdir/$outname.$common.*.bowtie2.hg19.$q.rc_insertion_points \
				| sort | uniq -c > $outdir/$outname.$common.both_ltr.bowtie2.hg19.$q.rc_insertion_points.rc_overlappers

		done	#	for q

		group_list=( )

	fi
done	#	for sample
