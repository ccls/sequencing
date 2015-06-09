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

	for fastabase in `find -L $1 -depth 1 -type f -name \*TCGA\*fasta | sed 's/\(^.*TCGA-[^-]*-[^-]*-\).*/\1/' | uniq` ; do
		echo $fastabase	#	/Volumes/box/working/output/HERV_K113_fasta/G49538.TCGA-14-1034-

		fasta_read_length_filter.sh -s r $fastabase*.pre_ltr.fasta
		fasta_read_length_filter.sh $fastabase*.post_ltr.fasta

		#	Hoping the reverse complements are on the proper side.

		herv_adjacent_comparison.sh $fastabase*.pre_ltr.trim*.fasta
		herv_adjacent_comparison.sh $fastabase*.post_ltr.trim*.fasta

#	base=${1%.*}		#	drop the extension
#	ext=${1##*.}		#	grab the extension
#	#	name=${base#*/}	#	just in case given path, drop the path
#
#	cmd="samtools view -S -h -F 4 $base.sam"
#	echo $cmd
#	$cmd | awk -v out=$base.tmp.counts '
#		( /^@SQ/ ){ ref[substr($2,4)] = 0 }
#		( !/^@/ ){ ref[$3]++ }
#		END{
#			for ( key in ref ) {
#				print key, ref[key] >> out
#			}
#		}'
#
#	#	awk options must be in specific order
#	#	(--posix needed to be after the -v's but is no longer needed anyway)
#
#	cat $base.tmp.counts | tr " " "\t" | tr "-" "\t" | sort -k2,2n > $base.counts
#	\rm $base.tmp.counts
#
#	other=$1

	done

	shift
done
