#!/usr/bin/env bash


function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` <query list file> <reference list file>"
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


#	$1 = query list
#	$2 = reference list

for line in `cat $1` ; do

	if [ `grep -c $line $2` -eq 0 ]; then
		echo $line
#		echo "not found"
#	else
#		echo "found"
	fi

done
