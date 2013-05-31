#!/bin/sh
#
#	split a fasta file into smaller fasta files with a 100000+ characters
#
#
#	would be nice to pass size
#
while [ $# -ne 0 ] ; do
	awk '
	BEGIN{
		file_number=0
		char_count=0
		f=sprintf("'$1'_%02d",++file_number)
	}
	{
		char_count+=length
		if(/^>/ && char_count >= 100000){
			f=sprintf("'$1'_%02d",++file_number)
			char_count=0
		}
		print>>f
	}' $1
	shift
done
