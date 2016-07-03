#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` --length INTEGER FOLDER(s)_CONTAINING_TCGA_FASTAS"
	echo
	echo "Example:"
	echo "  `basename $0` /Volumes/box/working/output/HERV_K113_fasta/"
	echo
	echo "Defaults:"
	echo "  --length/-l ..... : 50"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

length=50

while [ $# -ne 0 ] ; do
	case $1 in
		-l|--l*)
			shift; length=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done


while [ $# -ne 0 ] ; do
	echo $1

#	-depth will produce errors like the following on ec2. Use -maxdepth instead.
#	find: paths must precede expression: 2
#	Usage: find [-H] [-L] [-P] [-Olevel] [-D help|tree|search|stat|rates|opt|exec] [path...] [expression]
#	for fastabase in `find -L $1 -depth 1 -type f -name \*TCGA\*fasta | sed 's/\(^.*TCGA-[^-]*-[^-]*-\).*/\1/' | uniq` ; do
	for fastabase in `find -L $1 -maxdepth 1 -type f -name \*TCGA\*fasta | sed 's/\(^.*TCGA-[^-]*-[^-]*-\).*/\1/' | uniq` ; do
		echo $fastabase	#	/Volumes/box/working/output/HERV_K113_fasta/G49538.TCGA-14-1034-

		cmd="fasta_read_length_filter.sh -l $length -s r $fastabase*.pre_ltr.fasta"
		echo $cmd
		$cmd
		cmd="fasta_read_length_filter.sh -l $length $fastabase*.post_ltr.fasta"
		echo $cmd
		$cmd

		#	Hoping the reverse complements are on the proper side.

		cmd="herv_adjacent_comparison.sh $fastabase*.pre_ltr.trim$length.fasta"
		echo $cmd
		$cmd
		cmd="herv_adjacent_comparison.sh $fastabase*.post_ltr.trim$length.fasta"
		echo $cmd
		$cmd

	done

	shift
done
