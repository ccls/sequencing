#!/usr/bin/env bash
#
#	split fasta files into smaller fasta files with a 100000+ characters
#
#		fasta_size_split.sh my_file_1.fasta my_file_2.fasta
#
#	if first argument is a number, it will be used for the size
#


if [ $# -eq 0 ]; then
	echo
	echo "split fasta files into smaller fasta files with a 100000+ characters"
	echo
	echo "Usage:"
	echo
	echo "`basename $0` [optional max_characters_ish]"
	echo
	echo "Example:"
	echo "`basename $0` 10000 dna.fasta"
	echo
	exit
fi

tmp=`echo $1 | tr -cd '[:digit:]'`
# need the x's in case is blank
if [ "x${tmp}" == "x${1}" ] ; then
	size=$1
	shift
else
	size=100000
fi

while [ $# -ne 0 ] ; do
	base=${1%.*}	#	remove the extension

	now=`date "+%Y%m%d%H%M%S"`
	#	expecting trailing / later so make sure its here now...
	subdir=$base.${now}.pieces.nobackup/
	mkdir $subdir

	awk -v subdir="$subdir" -v base="$base" -v size="$size" '
	function reset(){
		char_count=0
		f=sprintf("%s%s.%06d.fasta",subdir,base,++file_number)
	}
	BEGIN{
		file_number=0
		reset()
	}
	{
		char_count+=length
		if(/^>/ && char_count >= size ){
			close(f)
			reset()
		}
		print>>f
	}
	END {
		close(f)
	}' $1
	shift
done
