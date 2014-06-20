#!/bin/sh

#
#	blastn_check.sh dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.20130725154032.pieces/*txt
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



if [ $# -eq 0 ]; then
	echo
	echo "checks the contents of the given blast files"
	echo
	echo "Usage:"
	echo
	echo "`basename $0` blast_output_file(s)"
	echo
	echo "Example:"
	echo "blastn_check.sh some_blastn*files.txt"
	echo
	exit
fi


#
#	It is unlikely that more than one will actually be called,
#	however, the while loop provides a nice wrapper.
#
while [ $# -ne 0 ] ; do
	blast_glob=$1

	#
	#	If there aren't any, this loop will crash
	#	ls: fake.fasta.20130802170906.pieces/*fasta_*.blastn.txt: No such file or directory
	#
	for blast in `ls $blast_glob` ; do
		echo "Checking $blast"
#
#	I think that double square brackets are only needed for regexps [[ ]]
#
		if [[ ! `head -1 $blast` =~ ^BLASTN ]] ; then
			echo "  *  First line in $blast is not correct"
		else
			echo 'First line is good (BLASTN)'
		fi
		if [[ ! `tail -1 $blast` =~ ^Gap\ Penalties ]] ; then
			echo "  *  Last line in $blast is not correct"
		else
			echo 'Last line is good (Gap Penalties)'
		fi

		echo 'Counting "BLASTN" lines'
		grep 'BLASTN' $blast | wc -l

		echo 'Counting "Gap Penalties" lines'
		grep 'Gap Penalties' $blast | wc -l

		echo 'Counting "Query= " lines'
		grep 'Query= ' $blast | wc -l

		echo 'Counting "Effective search space used:" lines'
		grep 'Effective search space used:' $blast | wc -l

#	If crashes, output is incomplete.  May not have "Effective search" line.  And "BLASTN" won't start on left.
#	Trying without anchor.


		echo "Checking for silent failures. Rerun if any found."
		echo "Checking for 'no longer exists' caused by database files temporarily disappearing?"
		grep -n "no longer exists in database" $blast

		echo "Checking for non-printable control characters (won't show the chars though) ..."
		grep -n '[[:cntrl:]]'  $blast

#		if [ `grep '^BLASTN' $blast | wc -l` -gt 1 ] ; then
#			echo "  *  Too many 'first lines' in $blast"
#		fi
#		if [ `grep '^Gap Penalties' $blast | wc -l` -gt 1 ] ; then
#			echo "  *  Too many 'last lines' in $blast"
#		fi
	done

	shift
done
