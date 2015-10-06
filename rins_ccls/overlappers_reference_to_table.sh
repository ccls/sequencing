#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo 
	echo "Usage:"
	echo
	echo "`basename $0` 'FIND PATTERN FOR OVERLAPPER_FILES'"
	echo
	echo "Note: execute from 1 directory above sample directory."
	echo
	echo "Example:"
	echo "`basename $0` '*ALL*overlappers_reference'"
	echo "`basename $0` '*Q20*overlappers_reference'"
	echo
	exit
fi

#TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers_reference

now=`date "+%Y%m%d%H%M%S"`




#while [ $# -ne 0 ] ; do
#	file=$1
#	attempting to deal with 
#		/bin/ls: Argument list too long.
for file in `find . -type f -depth 2 -name $1` ; do




#	echo
#	echo $file
	base=${file%.*}		#	drop the extension
#	echo $base
#	ext=${file##*.}		#	grab the extension	(overlappers or rc_overlappers)
#	echo $ext
	filename=${file##*/}	#	just in case given path
#	echo $filename
	subject=${filename%%.*}	#	delete everything after the first .  
# expecting filename like this ... TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
#	echo $subject

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	#[[ ${ext} =~ overlappers ]] && direction='F' || direction='R'	#	FAIL, both match
#	[[ ${ext} == overlappers ]] && direction='F' || direction='R'
#	echo $direction

#	TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
#	TCGA-41-5651-10A/TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.rc_insertion_points
#	rc_overlappers
#	TCGA-41-5651-10A.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.rc_insertion_points.rc_overlappers
#	TCGA-41-5651-10A

	#	what's --posix for? can't remember
	awk -v subject=$subject --posix '{
		printf "%s,%d,%s\n",$2,$1,subject
	}' $file >> tmpfile.$now

	shift
done


#chr5:80442266:R,228,TCGA-41-5651-10A
#chr7:107023897:R,1,TCGA-41-5651-10A
#chr7:2428619:R,2,TCGA-41-5651-10A
#chr7:2428996:R,1,TCGA-41-5651-10A
#chr9:42230337:R,16,TCGA-41-5651-10A


#	awk does not do multidimensional arrays
gawk -F, '{
	p[$1]++
	s[$3]++
	b[$1][$3]=$2
}
END{
	asorti(p)
	asorti(s)
	printf "position"
	for(subj in s)
		printf ",%s",s[subj]
	printf "\n"

	for(pos in p){
		printf p[pos]
		for(subj in s)
			printf ",%s",b[p[pos]][s[subj]]
		printf "\n"
	}
}' tmpfile.$now

#chr5:64388446:F TCGA-06-0125-10A 126
#chr5:64388446:F TCGA-06-0185-01A 48
#chr5:64388446:F TCGA-06-0211-10A 340
#chr5:64388446:F TCGA-06-0152-10A 128
#chr5:64388446:F TCGA-06-0185-10B 40

#		( /^@SQ/ ){ ref[substr($2,4)] = 0 }
#		( !/^@/ ){ ref[$3]++ }
#		END{
#			for ( key in ref ) {
#				print key, ref[key] >> out
#			}
#		}'
#


