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


dir=`dirname $0`
#
#	Why did I need gawk and not just awk?
#	I'm guessing its due to the character classes I used for matching.
#	In this gawk script, awk finds control characters on EVERY line
#	making it really kinda pointless.
#
gawk -v head="$head" -v tail="$tail" -f "$dir/blast_check.gawk" $@

