#!/usr/bin/env bash

#
#	pieces_check.sh  dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.20130725154032.pieces
#
#	count the contents of the given path
#
#	parse the blastn output files to ensure that they are complete and not truncated
#	or written to by multiple processes at the same time or anything like that.
#
#	a blast files should start with
#
#	BLASTN 2.2.25+
#	
#	
#	Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
#	Miller (2000), "A greedy algorithm for aligning DNA sequences", J
#	Comput Biol 2000; 7(1-2):203-14.
#	
#
#	... and end with ...
#
#
#  Database: Nucleotide collection (nt)
#    Posted date:  Mar 31, 2013  4:15 PM
#  Number of letters in database: 47,995,143,033
#  Number of sequences in database:  24,586,792
#
#
#
#	Matrix: blastn matrix 1 -2
#	Gap Penalties: Existence: 0, Extension: 2.5
#


script=`basename $0`

function usage(){
	echo
	echo "checks the contents of the given 'pieces' directory(ies)"
	echo
	echo "Usage:"
	echo
	echo "$script [--blast COMMAND] [--skip INTEGER] directory(ies)"
	echo
	echo "COMMANDS: blastn, tblastx"
	echo
	echo "Example:"
	echo "$script dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.20130725154032.pieces"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

blast='blastn'
skip=0
while [ $# -ne 0 ] ; do
	case $1 in
		-b|--b*)
			shift; blast=$1; shift ;;
		-s|--s*)
			shift;
			tmp=`echo $1 | tr -cd '[:digit:]'`
			if [ "x${tmp}" = "x${1}" ] ; then
				skip=$1
				shift
			else
				echo ; echo "skip value not an integer"
				usage
			fi ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

#
#	It is unlikely that more than one will actually be called,
#	however, the while loop provides a nice wrapper.
#
while [ $# -ne 0 ] ; do
	path=$1
	fasta_count=`ls $path/*.fasta | wc -l`
	blast_count=`ls $path/*.txt | wc -l`

	echo "Found $fasta_count fasta files"
	echo "Found $blast_count blast files"

	i=0

	#	20141126 - added -L to follow links (for nobackup on genepi)
	#	20150309 - removed fasta from expected name
	for blast_output in `find -L $path -type f -name \*.txt` ; do		#	change to 'find' loop?
		i=`expr $i + 1`
		if [ $i -le $skip ] ; then
			echo "skipping $blast_output"
			continue
		fi
		echo "Checking $blast_output"

		case $blast in
			'blastn')
				head="^BLASTN"
				tail="^Gap Penalties"
				;;
			'tblastx')
				head="^TBLASTX"
				tail="^Window for multiple hits"
				;;
			*)
				echo "Unrecognized blast command $blast"
		esac
		if [[ ! `head -1 $blast_output` =~ $head ]] ; then
			echo "  *  First line in $blast_output is not correct"
		else
			echo 'First line is good'
		fi
		if [[ ! `tail -1 $blast_output` =~ $tail ]] ; then
			echo "  *  Last line in $blast_output is not correct"
		else
			echo 'Last line is good'
		fi

		#	--blast doesn't make any difference in blast_check.sh ...... yet anyway
		blast_check.sh --blast $blast $blast_output

	done

	shift
done
