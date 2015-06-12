#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` bamfile1 bamfile2"
	echo
	echo "Example:"
	echo "  `basename $0` bamfile1 bamfile2"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ge 2 ] ; do
	echo $1
	echo $2

	base=${1%.*}		#	drop the extension
#	ext=${1##*.}		#	grab the extension
	name=${base#*/}	#	just in case given path

	cmd="samtools depth $1 $2"
	echo $cmd
	$cmd | awk '( $3!=0 && $4!=0 && ( ( $3>$4 && $3>(100*$4) ) || ( $4>$3 && $4>(100*$3) ) )){ print }'

	shift
	shift
done
