#!/bin/sh


script=`basename $0`

function usage(){
  echo
  echo "checks the contents of the given blast output file"
  echo
  echo "Usage:"
  echo
  echo "$script [--blast COMMAND] blast_output_file(s)"
  echo
  echo "COMMANDS: blastn, tblastx"
  echo
  echo "Example:"
  echo "$script dna/output/fallon_SFPB001A_filtered_20130722/trinity_input_single.fasta.blastn.txt"
  echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


blast='blastn'
while [ $# -ne 0 ] ; do
	case $1 in
		-b|--b*)
			shift; blast=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*)
			break;;
	esac
done

case $blast in
	'blastn' )
		head="BLASTN"
		tail="Gap Penalties"
		;;
	'tblastx' )
		head="TBLASTX"
		tail="Window for multiple hits"
		;;
	*)
		echo "Unrecognized blast command $blast"
esac


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

#	using 1 awk call rather than 6 greps.  MUCH FASTER

#	awk seems to filter out the control character problems that I'm searching for !!!!!!!!!
#	gawk (gnu's version of awk) WORKS!

#	I don't understand the requirement for the "" around $head and $tail.  Without them, they fail.
exec gawk '
(NR%100000 == 0){
	print "Read",FNR,"lines from",FILENAME >> "'$HOME/$script'.log"
}
(/'"$head"'/){ head++ }
(/'"$tail"'/){ tail++ }
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
	print "First line count :",head,":"
	print "Last line count :",tail,":"
	if( head != tail ){ print " * First and Lines lines are out of sync" }
	print "Query=  line count :",query,":"
	print "Effective search space used: line count :",effective,":"
	if( query != effective ){ print " * Query= and Effective search lines are out of sync" }
	print "no longer exists in database line count :",nolonger,":"
	if( nolonger > 0 ){ print " * TOO MANY no longer exists in database lines" }
	print "nonprint character line count :",nonprint,":"
	if( nonprint > 0 ){ print " * TOO MANY nonprintable characters" }
	print "control character line count :",control,":"
	if( control > 0 ){ print " * TOO MANY control characters" }
	print "---"
	close("'$HOME/$script'.log")
}' $@
