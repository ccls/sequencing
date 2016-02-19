#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo
	echo "No files given"
	echo
	echo "Usage:"
	echo
	echo "$0 BLASTN_OUTPUT.txt"
	echo
	echo "Example:"
	echo "$0 trinity_input_paired.blastn.txt"
	echo
	exit
fi

while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then

		dir=`dirname $0`
		awk -v num_alignments=2 -v num_descriptions=2 -f "$dir/trim_alignments_and_descriptions.awk" $1

	else
		echo "File '$1' not found?"
	fi
	shift
done


