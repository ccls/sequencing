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
	echo "Input bam or sam MUST be sorted by name."
	echo "Input bam or sam MUST contain mate info in flags."
	echo
	echo "Output is currently 2 fastq files."
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
	name=${base##*/}	#	just in case given path, drop the path

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

	#	read file
	#	buffer read if R1
	#	if R2 and seq name is same as buffered, print both
	#	MUST use gawk for the bitwise command "and"
	samtools view $flag $1 | gawk '
		( and( $2 , 64 ) ){
			b1=$1
			b10=$10
			b11=$11
		}
		( and( $2 , 128 ) ){
			if ( $1 == b1 ){
				if ( length(b10) != length($10) ){
					print $1 >> "'${base}.diff_length_reads'"
				}
				if ( length(b11) != length($11) ){
					print $1 >> "'${base}.diff_length_quality'"
				}
				print "@"b1"/1" >> "'${base}_R1.fastq'"
				print b10 >> "'${base}_R1.fastq'"
				print "+" >> "'${base}_R1.fastq'"
				print b11 >> "'${base}_R1.fastq'"

				print "@"$1"/2" >> "'${base}_R2.fastq'"
				print $10 >> "'${base}_R2.fastq'"
				print "+" >> "'${base}_R2.fastq'"
				print $11 >> "'${base}_R2.fastq'"
				b1=b10=b11=""
			}
		}'

	shift
done
