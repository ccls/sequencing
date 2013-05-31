#!/bin/sh
#
#	split fasta files into smaller fasta files with a 100000+ characters
#
#		fasta_size_split.sh my_file_1.fasta my_file_2.fasta
#
#	if first argument is a number, it will be used for the size
#

tmp=`echo $1 | tr -cd '[:digit:]'`
# need the x's in case is blank
if [ "x${tmp}" == "x${1}" ] ; then
	size=$1
	shift
else
	size=100000
fi

while [ $# -ne 0 ] ; do
	awk '
	BEGIN{
		file_number=0
		char_count=0
		f=sprintf("'$1'_%04d",++file_number)
	}
	{
		char_count+=length
		if(/^>/ && char_count >= '$size'){
			close(f)
			f=sprintf("'$1'_%04d",++file_number)
			char_count=0
		}
		print>>f
	}' $1
	shift
done
