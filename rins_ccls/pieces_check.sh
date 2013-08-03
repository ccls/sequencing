#!/bin/sh

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

while [ $# -ne 0 ] ; do
	path=$1

#	re-adding .fasta extension to increase findability.
#	also may change the number of digits in pieces.
#	would prefer to be able to glob like a regex \d+ , but can't
#
#	fasta_count=`ls $path/*fasta_????????? | wc -l`
#	fake.fasta_000000006.fasta
	fasta_count=`ls $path/*fasta_*.fasta | wc -l`
#	blast_count=`ls $path/*fasta_?????????.blastn.txt | wc -l`
#	fake.fasta_000000006.fasta.blastn.txt
	blast_count=`ls $path/*fasta_*.blastn.txt | wc -l`

	echo "Found $fasta_count fasta files"
	echo "Found $blast_count blast files"

#	for blast in `ls $path/*fasta_?????????.blastn.txt` ; do
	for blast in `ls $path/*fasta_*.blastn.txt` ; do
		echo "Checking $blast"
#
#	I think that double square brackets are only needed for regexps [[ ]]
#
		if [[ ! `head -1 $blast` =~ ^BLASTN ]] ; then
			echo "  *  First line in $blast is not correct"
#		else
#			echo 'good head'
		fi
		if [[ ! `tail -1 $blast` =~ ^Gap\ Penalties ]] ; then
			echo "  *  Last line in $blast is not correct"
#		else
#			echo 'good tail'
		fi

		if [ `grep '^BLASTN' $blast | wc -l` -gt 1 ] ; then
			echo "  *  Too many 'first lines' in $blast"
		fi
		if [ `grep '^Gap Penalties' $blast | wc -l` -gt 1 ] ; then
			echo "  *  Too many 'last lines' in $blast"
		fi
	done

	shift
done


