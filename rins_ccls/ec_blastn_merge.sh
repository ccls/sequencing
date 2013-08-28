#!/bin/sh

#	ec_blastn_merge.sh *.pieces

#	each filename on the command line
while [ $# -ne 0 ] ; do
	if [ -d $1 ] ; then

		echo "concatenating $1/*blastn.txt to $1.blastn.txt"

		if [ -f $1.blastn.txt ] ; then
			echo "$1.blastn.txt already exists.  Skipping."
		else
			cat $1/*blastn.txt > $1.blastn.txt
		fi

	else
		echo "$1 doesn't seem to be a directory."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
