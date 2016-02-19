#!/usr/bin/env bash
#
# split fasta files into smaller fasta files with no more than 200 reads
#
#		fasta_read_split.sh my_file_1.fasta my_file_2.fasta
#
#	if first argument is a number, it will be used for the max_reads
#
#	
if [ $# -eq 0 ]; then
	echo 
	echo "split fasta files into smaller fasta files with no more than 200 reads"
	echo
	echo "Usage:"
	echo
	echo "`basename $0` [optional max_read_size]"
	echo
	echo "Example:"
	echo "`basename $0` 100 dna.fasta"
	echo
	exit
fi

tmp=`echo $1 | tr -cd '[:digit:]'`
#	need the x's in case is blank
if [ "x${tmp}" == "x${1}" ] ; then
	max_reads=$1
	shift
else
	max_reads=200
fi

while [ $# -ne 0 ] ; do
	base=${1%.*}	#	remove the extension

	now=`date "+%Y%m%d%H%M%S"`
	#	expecting trailing / later so make sure its here now...
	subdir=$base.${now}.pieces.nobackup/
	mkdir $subdir

	awk -v subdir="$subdir" -v base="$base" -v max_reads="$max_reads" '
	function reset(){
		read_count=0
		f=sprintf("%s%s.%06d.fasta",subdir,base,++file_number)
	}
	BEGIN{
		file_number=0
		reset()
	}	
	{
		if(/^>/){
			if( read_count >= max_reads ){
				close(f)
				reset()
			}
			read_count++
		}
		print>>f
	}
	END {
		close(f)
	}' $1
	shift
done
