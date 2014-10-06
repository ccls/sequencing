#!/bin/sh

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
#	moved comments outside of awk program so don't show in 'ps -ef'




		awk '
		{
			if( $1 == "Query=" ){
				align_count=0
				desc_count=-1
			}
			if( /^>/ ){
				align_count++
				desc_count=-1
			}
			if( /^Lambda     K      H/ ){
				align_count=0
			}


			if( desc_count >= 0 ){
				if( !/^\s*$/ ){
					desc_count++
				}
			}
			if( desc_count > 3 ){
				if( /^\s*$/ ){
					desc_count=-1
				}
			}
			if( /^Sequences producing significant alignments/ ){
				desc_count=0
			}



			if( ( align_count <= 2 ) && ( desc_count <= 2 ) ){
				print $0
			}
		}' $1

	else
		echo "File '$1' not found?"
	fi
	shift
done


