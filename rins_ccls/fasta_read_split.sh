#!/bin/sh
#
# split fasta files into smaller fasta files with no more than 200 reads
#
#		fasta_read_split.sh my_file_1.fasta my_file_2.fasta
#
#	if first argument is a number, it will be used for the max_reads
#

tmp=`echo $1 | tr -cd '[:digit:]'`
#	need the x's in case is blank
if [ "x${tmp}" == "x${1}" ] ; then
	max_reads=$1
	shift
else
	max_reads=200
fi

while [ $# -ne 0 ] ; do
	awk '
	BEGIN{
		file_number=0
		read_count=0
		f=sprintf("'$1'_%04d",++file_number)
	}
	{
		if(/^>/){
			if( read_count >= '$max_reads' ){
				close(f)
				f=sprintf("'$1'_%04d",++file_number)
				read_count=0
			}
			read_count++
		}
		print>>f
	}' $1
	shift
done
