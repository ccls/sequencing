#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` tumor.fasta normal.fasta"
	echo
	echo "Example:"
	echo "  `basename $0` tumor.fasta normal.fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -ne 2 ] && usage


echo $1, $2
base=${1%.*}		#	drop the extension
ext=${1##*.}		#	grab the extension
#	name=${base#*/}	#	just in case given path, drop the path

fastx_collapser -i $1 -o $base.collapsed.$ext

bowtie2-build $base.collapsed.$ext $base

bowtie2 -f -x $base -U $2 -S $base.sam
\rm *bt2

cmd="samtools view -S -h -F 4 $base.sam"
echo $cmd
$cmd | awk -v out=$base.counts '
	( !/^@/ ){ ref[$3]++ }
	END{
		for ( key in ref ) {
			print key, ref[key] >> out
		}
	}'

#	awk options must be in specific order 
#	(--posix needed to be after the -v's but is no longer needed anyway)
