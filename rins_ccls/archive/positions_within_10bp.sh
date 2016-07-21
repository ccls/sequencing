#!/usr/bin/env bash


function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` file1 file2"
	echo
	echo "Expecting file content in the format of ..."
	echo
	echo "..."
	echo "chr8:99063786"
	echo "chr9:105365407"
	echo "chrX:74554211"
	echo "chrY:14564844"
	echo "..."
	echo
	echo "Example:"
	echo "  `basename $0` file1 file2"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -ne 2 ] && usage
if [ ! -f $1 ] ; then 
	echo "$1 is not a file" 
	usage
fi
if [ ! -f $2 ] ; then
	echo "$2 is not a file" 
	usage
fi


#	$1 = queries	
#	$2 hg19 reference

for line in `cat $1` ; do

#	echo $line
	chr=${line%%:*}	#	remove everything after first colon (including colon)
	line=${line#*:}	#	remove everything before first colon (including colon)
	pos=${line%%:*}	#	remove everything after first colon (including colon)
#	echo $chr
#	echo $pos

	#	Expecting file content format like so ...
	#	chrY:6616930:

	awk -F: -v chr="$chr" -v pos="$pos" '
		( ( $1 == chr ) && ( (pos-10) < $2 ) && ( (pos+10) > $2 ) ){
			print
		}' $2

done

