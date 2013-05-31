#!/bin/sh
#
# split a fasta file into smaller fasta files with no more than 200 reads
#
#
#	would be nice to pass a read number
#
while [ $# -ne 0 ] ; do
	awk '
	BEGIN{
		file_number=0
		read_count=0
		f=sprintf("'$1'_%02d",++file_number)
	}
	{
		if(/^>/){
			if( read_count >= 200 ){
				f=sprintf("'$1'_%02d",++file_number)
				read_count=0
			}
			read_count++
		}
		print>>f
	}' $1
	shift
done
