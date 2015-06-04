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
other=$2
while [ $# -ne 0 ] ; do

	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
	#	name=${base#*/}	#	just in case given path, drop the path

	cmd="fastx_collapser -i $1 -o $base.collapsed.$ext"
	echo $cmd
	$cmd

	cmd="bowtie2-build $base.collapsed.$ext $base"
	echo $cmd
	$cmd

	cmd="bowtie2 -f -x $base -U $other -S $base.sam"
	echo $cmd
	$cmd
	\rm $base*bt2

	#	-F 4 will stop the *'s from coming out at night.

	cmd="samtools view -S -h -F 4 $base.sam"
	echo $cmd
	$cmd | awk -v out=$base.tmp.counts '
		( /^@SQ/ ){ ref[substr($2,4)] = 0 }
		( !/^@/ ){ ref[$3]++ }
		END{
			for ( key in ref ) {
				print key, ref[key] >> out
			}
		}'

	#	awk options must be in specific order
	#	(--posix needed to be after the -v's but is no longer needed anyway)

	cat $base.tmp.counts | tr " " "\t" | tr "-" "\t" | sort -k2,2n > $base.counts
	\rm $base.tmp.counts

	other=$1
	shift
done
