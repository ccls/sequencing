#!/bin/sh

function usage(){
	echo
	echo "Creates a directory of the name of each of the given fasta files and appends"
	echo "the date and .pieces.  The fasta file is then split into fasta files each"
	echo "with no more than the maximum number of reads and sequentially numbers them"
	echo "placing them in the aforementioned directories."
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--command COMMAND] [--max_reads INTEGER] [--dbs COMMA_SEP_STRING] [--outfmt BLASTN_OUTFMT#] fasta_filelist"
	echo
	echo "The default command is 'blastn'."
	echo
	echo "The default max reads per piece is 1000."
	echo
	echo "The default dbs are just nt."
	echo 
	echo "The default outfmt is 0."
	echo 
	echo "Example: `basename $0` -m 500 --dbs nt,viral,hg /my/path/*fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[[ $# -eq 0 ]] && usage

uname=`uname -n`
command='blastn'
dbs='nt'

#       leading with the ": " stops execution
#       just ${BLASTDB:"/Volumes/cube/working/indexes"}
#       would try to execute the result.  I just want the OR/EQUALS feature
: ${BLASTDB:="/Volumes/cube/working/indexes"}

#       they MUST be exported, apparently, to be picked up by calls
export BLASTDB

#if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
#	BLASTDB='/my/home/jwendt/dna/blast/nt'
#fi

max_reads=1000
outfmt=0
while [ $# -ne 0 ] ; do
	case $1 in
		-c|--c*)
			shift; command=$1; shift ;;
		-d|--d*)
			shift; dbs=$1; shift ;;
		-m|--m*)
			shift; 
			tmp=`echo $1 | tr -cd '[:digit:]'`
			if [ "x${tmp}" = "x${1}" ] ; then
				max_reads=$1
				shift
			else
				#	if max reads isn't an integer, awk works oddly
				echo ; echo "max reads value not an integer"
				usage
			fi ;;
		-ou|--ou*)
			shift; outfmt=$1; shift ;;
		-op|--op*)
			shift; options=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done

#	each filename on the command line
while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then
		now=`date "+%Y%m%d%H%M%S"`

		fasta_base=`basename $1 | sed 's/\(.*\)\..*/\1/'`

		#	extension=${filename##*.}
		#	the longest matching  pattern  (the  ``##'' case)  deleted
		#	It deletes the longest *. match leaving only the extension behind (starts at the beginning (left))
		#	% and %% start from the trailing or right side (could be similar to using basename)

		#	not used so why bother
		#		fasta_dir=`dirname $1`

		subdir=`basename $1`.${now}.pieces

		#
		#	if given path to file, put output there or where command was executed???
		#	I say where command is executed.
		#
		mkdir $subdir

		#
		#	I would like to drop the original fasta file's extension,
		#	but that would require changes below.  Be careful.
		#	Keeping the .fasta in the blast name does help clarify which have run
		#	and which have not, so I say keep if for the moment.
		#

		#
		#	Changed from 8 digits to 6.  Only ever had a 15000 max, so even 5 would be realistic.
		#	Perhaps count the reads in the file like I do in rins_ccls/bioruby_lane_fasta.rb
		#	if need to be more specific.
		#		max_digits = Math.log10( total_sequences + 1 ).ceil
		#

		#
		#	I'd also like to undo the duplication of setting the read_count, f and file_number
		#

		awk '
			BEGIN{
				file_number=0
				read_count=0
				f=sprintf("'$subdir/$fasta_base'.%06d.fasta",++file_number)
			}
			{
				if(/^>/){
					if( read_count >= '$max_reads' ){
						close(f)
						f=sprintf("'$subdir/$fasta_base'.%06d.fasta",++file_number)
						read_count=0
					}
					read_count++
				}
				print>>f
			}' $1

		negative_gilist=''
		if [ -f 'negative_gilist' ] ; then
			negative_gilist="-negative_gilist $PWD/negative_gilist"
		fi

		#	on some occassions, this list is too long for ls so changing to find
		#	for file in `ls $PWD/$subdir/${fasta_base}_*.fasta` ; do
		#	interesting _*. is ok, but .*. is not.  must escape the * here so .\*.
		for file in `find $PWD/$subdir -type f -name ${fasta_base}.\*.fasta` ; do

			#	allowing for multiple db blasting
			for db in `echo $dbs | sed 's/,/ /g'` ; do

				db_base_name=`basename $db`

				#
				# BE ADVISED!  For whatever reason, the order of some options does matter.
				#              I don't know why, but if negative_gilist is last, it is ignored.
				#              Could be others.
				#

				#cmd="$command $negative_gilist -show_gis -query $file -db $db \
				cmd="$command $negative_gilist -query $file -db $db \
					-num_alignments 20 -num_descriptions 30 -evalue 0.05 -outfmt $outfmt \
					-out $file.blastn_${db_base_name}.txt $options"

				#	20140724 - added -show_gis to potentially help with this Uncultured stuff
				#	20140729 - added -num_descriptions 30 to help minimize file size
				#	20140825 - commented out show_gis as may be contributing to larger output size

				echo $cmd

				#	20140729 - removed eval wrapper around the push command
				#eval "simple_queue.sh push '$cmd'"
				simple_queue.sh push "$cmd"

			done

		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
