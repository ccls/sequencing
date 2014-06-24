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

#			echo 'Counting "BLASTN" lines'
#			grep 'BLASTN' $blast | wc -l
#	
#			echo 'Counting "Gap Penalties" lines'
#			grep 'Gap Penalties' $blast | wc -l
#	
#			echo 'Counting "Query= " lines'
#			grep 'Query= ' $blast | wc -l
#	
#			echo 'Counting "Effective search space used:" lines'
#			grep 'Effective search space used:' $blast | wc -l
#	
#	#	If crashes, output is incomplete.  May not have "Effective search" line.  And "BLASTN" won't start on left.
#	#	Trying without anchor.
#	
#	
#			echo "Checking for silent failures. Rerun if any found."
#			echo "Checking for 'no longer exists' caused by database files temporarily disappearing?"
#			grep -n "no longer exists in database" $blast
#	
#	
#	#	not sure what the diffs would be here for these POSIX char classes
#	#    [[:ascii:]] - matches a single ASCII char
#	#    [^[:ascii:]] - matches a single non-ASCII char
#	#[^[:print:]] will probably suffice for you.**
#	
#			echo "Checking for non-printable control characters (won't show the chars though) ..."
#			grep -n '[[:cntrl:]]'  $blast
#	
#	#		if [ `grep '^BLASTN' $blast | wc -l` -gt 1 ] ; then
#	#			echo "  *  Too many 'first lines' in $blast"
#	#		fi
#	#		if [ `grep '^Gap Penalties' $blast | wc -l` -gt 1 ] ; then
#	#			echo "  *  Too many 'last lines' in $blast"
#	#		fi







# Here are the character classes defined by the POSIX standard.
#
#[:alnum:]
#    Alphanumeric characters. 
#[:alpha:]
#    Alphabetic characters. 
#[:blank:]
#    Space and tab characters. 
#[:cntrl:]
#    Control characters. 
#[:digit:]
#    Numeric characters. 
#[:graph:]
#    Characters that are printable and are also visible. (A space is printable, but not visible, while an `a' is both.) 
#[:lower:]
#    Lower-case alphabetic characters. 
#[:print:]
#    Printable characters (characters that are not control characters.) 
#[:punct:]
#    Punctuation characters (characters that are not letter, digits, control characters, or space characters). 
#[:space:]
#    Space characters (such as space, tab, and formfeed, to name a few). 
#[:upper:]
#    Upper-case alphabetic characters. 
#[:xdigit:]
#    Characters that are hexadecimal digits. 


#	interesting note.  [:cntrl] in awk includes carriage returns (not so in grep).  using non-printable

#	http://sed.sourceforge.net/sedfaq3.html
#     [[:cntrl:]]  - [\x00-\x19\x7F] Control characters
#			(/[:cntrl:][^]/){ 

		#	using 1 awk call rather than 6 greps.  MUCH FASTER

		#	awk seems to filter out the control character problems that I'm searching for !!!!!!!!!
		#	gawk (gnu's version of awk) WORKS!

		gawk '
			(/BLASTN/){ blastn++ }
			(/Gap Penalties/){ gap++ }
			(/Query= /){ query++ }
			(/Effective search space used:/){ effective++ }
			(/no longer exists in database/){ nolonger++ }
			(/[^[:print:]]/){ 	#	too many (includes \anything?)
				print "Non-printable character(s) found on line number :",NR,":"
				print $0
				nonprint++ 
			}
			(/[[:cntrl:]]/){	#	REQUIRES double brackets as they are for different reasons
				print "Control character(s) found on line number :",NR,":"
				print $0
				control++ 
			}
			END{
				print "BLASTN line count"
				print blastn
				print "Gap Penalties line count"
				print gap
				if( blastn != gap ){ print " * BLASTN and Gap Penalties lines are out of sync" }
				print "Query=  line count"
				print query
				print "Effective search space used: line count"
				print effective
				if( query != effective ){ print " * Query= and Effective search lines are out of sync" }
				print "no longer exists in database line count"
				print nolonger
				if( nolonger > 0 ){ print " * TOO MANY" }
				print "nonprint character line count"
				print nonprint
				if( nonprint > 0 ){ print " * TOO MANY" }
				print "control character line count"
				print control
				if( control > 0 ){ print " * TOO MANY" }
				print "---"
			}
		' $blast

	done	#	for blast in `ls $blast_glob` ; do

	shift
done
