#!/usr/bin/env bash
#	manually set -x as can't pass through env
set -x

###!/bin/sh -x

#
#	pieces_merge.sh  dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.20130725154032.pieces
#


script=`basename $0`

function usage(){
	echo
	echo "merge the blast contents of the given 'pieces' directory(ies)"
	echo
	echo "Usage:"
	echo
	echo "$script directory(ies)"
	echo
	echo "Example:"
	echo "$script dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.20130725154032.pieces"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

touch md5sums
chmod +w md5sums

#
#	It is unlikely that more than one will actually be called,
#	however, the while loop provides a nice wrapper.
#
while [ $# -ne 0 ] ; do

	#	It is possible that a single directory is used for multiple blasting.
	#	Unlikely.  But possible.
	#	trinity_non_human_single.fasta.20150125184142.pieces.nobackup/trinity_non_human_single.000001.fasta.blastn_viral_genomic.txt
	#	XXX is just used as a placeholder.
	for pieces in `find -L $1 -type f -name \*.txt | sed 's/\.[0-9]\{6\}\./.XXX./g' | uniq` ; do

		date
		echo $pieces
		#	trinity_non_human_single.fasta.20150125184145.pieces.nobackup/trinity_non_human_single.XXX.blastn_nt.txt
		base=`basename $pieces`
		#	trinity_non_human_single.XXX.blastn_nt.txt

		outname=${base/.XXX./.}
		outname=${outname/.fasta./.}

		mmv -a "${pieces/.XXX./*}" $outname
		#mmv -a "trinity_non_human_single.fasta.20150125184145.pieces.nobackup/trinity_non_human_single*blastn_nt.txt" trinity_non_human_single.blastn_nt.txt

		chmod -w $outname
		md5sum $outname >> md5sums
		gzip --best  $outname
		md5sum $outname.gz >> md5sums

	done
	shift
done

chmod -w md5sums

