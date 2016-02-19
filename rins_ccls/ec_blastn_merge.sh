#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo
	echo "Loops through given directories, concatenating the contents of each blastn.txt"
	echo "into a blastn.txt of the same name as each given directory."
	echo
	echo "Usage:"
	echo
	echo "`basename $0` directory(ies)"
	echo
	echo "Example: `basename $0` *.pieces"
	echo
	exit	#	not really necessary as loop won't happen
fi

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
