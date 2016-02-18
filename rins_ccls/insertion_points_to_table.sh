#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo 
	echo "Usage:"
	echo
	echo "`basename $0` 'FIND PATTERN FOR INSERTION POINT FILES'"
	echo
	echo "Example:"
	echo "`basename $0` '*Q00*points'"
	echo "`basename $0` '*Q20*points'"
	echo
	exit
fi

#TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.rc_insertion_points

now=`date "+%Y%m%d%H%M%S"`

tmpfile="tmpfile.$1.$now"


#while [ $# -ne 0 ] ; do
#	file=$1
#	attempting to deal with 
#		/bin/ls: Argument list too long.
#	Use find instead of expanding argument list (so will need quoted)
for file in `find ./ -type f -depth 2 -name $1` ; do

#	echo
#	echo $file
	base=${file%.*}		#	drop the extension
#	echo $base
	ext=${file##*.}		#	grab the extension	(insertion_points or rc_insertion_points)
#	echo $ext
	filename=${file##*/}	#	just in case given path
#	echo $filename
	subject=${filename%%.*}	#	delete everything after the first .  
# expecting filename like this ... TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.rc_insertion_points
#	echo $subject

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	#[[ ${ext} =~ insertion_points ]] && direction='F' || direction='R'	#	FAIL, both match
	[[ ${ext} == insertion_points ]] && direction='F' || direction='R'
#	echo $direction

#	TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
#	TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.rc_insertion_points
#	rc_overlappers
#	TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
#	TCGA-41-5651-10A


	#	There ALSO will be pre and post for insertion points
	#		*.post_ltr.*.insertion_points
	#		*.post_ltr.*.rc_insertion_points
	#		*.pre_ltr.*.insertion_points
	#		*.pre_ltr.*.rc_insertion_points

	[[ ${filename} =~ pre_ltr ]] && pre_or_post='PRE' || pre_or_post='POST'

	#	Expecting file content like ...
	#	chr4:179554141
	#	chr4:179554141
	#	chr4:179554141
	#	chr4:179554141
	#	chr4:208766
	#	chr4:25249260
	#	chr4:3094365
	#	chr4:3094365

	#	uniq -c will convert the above to ...
	#		4	chr4:179554141
	#		1	chr4:208766
	#		1	chr4:25249260
	#		2	chr4:3094365

	#	uniq -c $file | awk -v subject=$subject -v direction=$direction \
	#	The input file needs to be grouped.  Should already be sorted, but in case.
	#	Sort will sort differently but still grouped so OK.
	sort $file | uniq -c | awk -v subject=$subject -v direction=$direction \
		-v pre_or_post=$pre_or_post --posix '{
		printf "%s:%s:%s,%d,%s\n",$2,direction,pre_or_post,$1,subject
	}' >> $tmpfile

	shift
done


#chr5:80442266:R,228,TCGA-41-5651-10A
#chr7:107023897:R,1,TCGA-41-5651-10A
#chr7:2428619:R,2,TCGA-41-5651-10A
#chr7:2428996:R,1,TCGA-41-5651-10A
#chr9:42230337:R,16,TCGA-41-5651-10A


dir=`dirname $0`
gawk -f "$dir/to_table.gawk" $tmpfile

