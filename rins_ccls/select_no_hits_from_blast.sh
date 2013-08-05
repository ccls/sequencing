#!/bin/sh

#	Select the sequence name / line from the blast output file
#	where there are "No hits found"

#
#	This would be perfect IF it was always 6. 
#	The sequence name can span multiple lines mostly due to complex paths..
#	It is usually 1 or 2, but not always.
#
#grep --before-context=6 'No hits found' trinity_non_human_paired.blastn.txt | head


if [ $# -eq 0 ]; then
	echo
	echo "No files given"
	echo
	echo "Usage:"
	echo
	echo "select_no_hits_from_blast.sh BLASTN_OUTPUT.txt > LIST_OF_NO_HITS_FOUND.names &"
	echo
	echo "Example:"
	echo "select_no_hits_from_blast.sh trinity_input_paired.blastn.txt > trinity_input_paired.blastn.no_hits_found.names &"
	echo
fi
while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then
		awk '
		BEGIN{
			seq_name=""
			seq_line=""
		}
		{
#
#	Query= comp102_c0_seq1 len=104 path=[1:0-82 84:83-103]
#
#	if(/^Query= \w+ /){
#
#	... works for above, but not with ...
#
#	Query= HS3:309:D2385ACXX:2:1101:1249:2272/1
#
#			if(/^Query= \S+ /){
#	that does not seem to work either, but this does 
#	( as awk breaks line into columns on whitespace )
#
			if( $1 == "Query=" ){
				seq_name=$2
				seq_line=$0
			} else if(/No hits found/){
				print seq_name
#				print seq_line
			}
		}' $1
	else
		echo "File '$1' not found?"
	fi
	shift
done

