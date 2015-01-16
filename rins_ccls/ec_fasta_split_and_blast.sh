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
	echo "`basename $0` [--command COMMAND] [--max_reads INTEGER] [--dbs COMMA_SEP_STRING] "
	echo "              [--evalue FLOAT] [--dust STRING]"
	echo "              [--num_alignments INTEGER] [--num_descriptions INTEGER]"
	echo "              [--prefix STRING] [--suffix STRING]"
	echo "              [--std_out_only]"
	echo "              [--outfmt BLASTN_OUTFMT#] fasta_filelist"
	echo
	echo "Notes:"
	echo "  The dust value will require quotes if contains spaces like the default."
	echo "  (blastx and tblastx do not use -dust)"
	echo
	echo "Defaults:"
	echo "  command . : blastn"
	echo "  max_reads : 1000"
	echo "  dbs ..... : nt"
	echo "  outfmt .. : 0"
	echo "  evalue .. : 0.05"
	echo "  dust .... : ''	(blastn default is '20 64 1')"
	echo "  num_alignments .... : 20"
	echo "  num_descriptions .. : 20"
	echo "  prefix .. : "
	echo "  suffix .. : "
	echo 
	echo "Example:"
	echo "  `basename $0` --max 500 --dbs nt,viral,hg /my/path/*fasta"
	echo "  `basename $0` --dust '10 32 1' /my/path/*fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

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
evalue=0.05
num_alignments=20
num_descriptions=20
prefix=''
suffix=''
std_out_only='false'
dust_value=""	#	"20 64 1"
while [ $# -ne 0 ] ; do
	case $1 in
		-c|--c*)
			shift; command=$1; shift ;;
		--db*)
			shift; dbs=$1; shift ;;
		--e*)
			shift; evalue=$1; shift ;;
		--du*)
			shift; dust_value=$1; shift ;;
		--num_a*)
			shift; num_alignments=$1; shift ;;
		--num_d*)
			shift; num_descriptions=$1; shift ;;
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
		--pr*)
			shift; prefix=$1; shift ;;
		--su*)
			shift; suffix=$1; shift ;;
		--st*)
			std_out_only='true'; shift ;;
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
		#	20150116 - using
		fasta_dir=`dirname $1`

		#subdir=`basename $1`.${now}.pieces
		subdir=$fasta_dir/`basename $1`.${now}.pieces

		#
		#	if given path to file, put output there or where command was executed???
		#	I say where command is executed.
		#	20150116 - And now I say where the file is.
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

		dust=''
		if [ -n "$dust_value" ] ; then	#	quotes are required
			dust="-dust ''$dust_value''"
		fi

		#	on some occassions, this list is too long for ls so changing to find
		#	for file in `ls $PWD/$subdir/${fasta_base}_*.fasta` ; do
		#	interesting _*. is ok, but .*. is not.  must escape the * here so .\*.
#		for file in `find $PWD/$subdir -type f -name ${fasta_base}.\*.fasta` ; do
		for file in `find $subdir -type f -name ${fasta_base}.\*.fasta` ; do

			#	allowing for multiple db blasting
			for db in `echo $dbs | sed 's/,/ /g'` ; do

				db_base_name=`basename $db`

				#
				# BE ADVISED!  For whatever reason, the order of some options DOES matter.
				#              I don't know why, but if negative_gilist is last, it is ignored.
				#              Could be others.
				#

				cmd="$prefix $command $negative_gilist -query $file -db $db $dust \
					-num_alignments $num_alignments -num_descriptions $num_descriptions \
					-evalue $evalue -outfmt $outfmt \
					-out $file.${command}_${db_base_name}.txt $options $suffix"


				#	20140724 - Added -show_gis to potentially help with this Uncultured stuff.
				#	20140729 - Added -num_descriptions 30 to help minimize file size.
				#	20140825 - Commented out show_gis as may be contributing to larger output size.
				#	20141003 - Added num_alignments, num_descriptions, evalue and dust.
				#            The dust value needs double single quoted when simple queue pushes it.

				echo $cmd

				if [ $std_out_only != 'true' ] ; then
					#	20140729 - removed eval wrapper around the push command
					#eval "simple_queue.sh push '$cmd'"
					simple_queue.sh push "$cmd"
				fi

			done

		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
