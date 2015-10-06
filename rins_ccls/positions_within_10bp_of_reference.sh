#!/usr/bin/env bash


function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` referencefile samplefilelist"
	echo
	echo "Expecting file content in the format of ..."
	echo
	echo "..."
	echo "chr8:99063786:F"
	echo "chr9:105365407:F"
	echo "chrX:74554211:R"
	echo "chrY:14564844:R"
	echo "..."
	echo
	echo "With or without the direction is acceptable."
	echo
	echo "Example:"
	echo "  `basename $0` referencefile samplefilelist"
	echo
	echo "find . -type d -depth 1 -execdir sh -c 'echo {}; cd {}; positions_within_10bp_of_reference.sh ../../overlapper_reference.Q20 *Q20*ts | sort | uniq -c > {}.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.hg19.Q20.overlappers_reference' \;"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
#[ $# -ne 2 ] && usage
[ $# -lt 2 ] && usage
if [ ! -f $1 ] ; then 
	echo "$1 is not a file" 
	usage
fi
if [ ! -f $2 ] ; then
	echo "$2 is not a file" 
	usage
fi

reference=$1
shift

while [ $# -ne 0 ] ; do
#	echo $1
	for line in `cat $reference` ; do

		#	CORRECTION .... chrY:4395088:F ... REFERENCE INCLUDES DIRECTION

#		echo $line
		chr=${line%%:*}	#	remove everything after first colon (including colon)
		line2=${line#*:}	#	remove everything before first colon (including colon)
		pos=${line2%%:*}	#	remove everything after first colon (including colon)
#		echo $chr
#		echo $pos

		#	Expecting file content format like so ...
		#	chrY:6616930:

#	print the reference insertion point or the point found?

		awk -F: -v chr="$chr" -v pos="$pos" -v line="$line" '
			( ( $1 == chr ) && ( (pos-10) < $2 ) && ( (pos+10) > $2 ) ){
				print line
			}' $1
#				print line" - "$0

	done	#	for line in reference
	shift	#	change sample input file
done	#	while more input files
