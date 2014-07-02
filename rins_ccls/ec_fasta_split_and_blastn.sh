#!/bin/sh

function usage(){
	echo
	echo "Creates a directory of the name of each of the given fasta files and appends"
	echo "the date and .pieces.  The fasta file is then split into fasta files each"
	echo "with no more than the maximum number of reads and sequentially numbers them"
	echo "placing them in the aforementioned directories."
	echo
	echo "Usage:"
	echo
	echo "`basename $0` [--max_reads INTEGER] [--dbs COMMA_SEP_STRING] fasta_filelist"
	echo
	echo "The default max reads per piece is 1000."
	echo
	echo "The default dbs are just nt."
	echo 
	echo "Example: `basename $0` -m 500 --dbs nt,viral,hg /my/path/*fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[[ $# -eq 0 ]] && usage

uname=`uname -n`

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
while [ $# -ne 0 ] ; do
	case $1 in
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
		-o|--o*)
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


		#	on some occassions, this list is too long for ls so changing to find
		#	for file in `ls $PWD/$subdir/${fasta_base}_*.fasta` ; do
		#	interesting _*. is ok, but .*. is not.  must escape the * here so .\*.
		for file in `find $PWD/$subdir -type f -name ${fasta_base}.\*.fasta` ; do

			#	allowing for multiple db blasting
			for db in `echo $dbs | sed 's/,/ /g'` ; do
				#
				#	No longer putting cluster stuff in the queued command.
				#	This way, I can use it whereever I like.
				#
#				cmd=''	#	gotta reset it
#				if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
#					#	trinity_input_single.uniq.fasta_00000416.fasta
#					#	=> 'uniq' for num.  oops
#					#	num=`basename $file | awk -F. '{print $2}' | awk -F_ '{print $NF}'`
#					#	num=`basename $file | awk -F. '{print $(NF-1)}' | awk -F_ '{print $NF}'`
#					num=`basename $file | awk -F. '{print $(NF-1)}'`
#					cmd="srun --share --job-name=${num}_$db"
#				fi
				#echo db $db
#				cmd="$cmd blastn -query $file -db $db -num_alignments 20 -evalue 0.05 -outfmt 0 -out $file.blastn_${db}.txt $options &"

				db_base_name=`basename $db`
				cmd="blastn -query $file -db $db -num_alignments 20 -evalue 0.05 -outfmt 0 -out $file.blastn_${db_base_name}.txt $options"
				#cmd="$cmd blastn_wrapper.sh $file $db &"

				echo $cmd

#				if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
					#	need to eval to use the &
					#	want the & in the queue'd command, not here.
					#	if were here, will cause database error by trying to write to it at same time
					eval "simple_queue.sh push '$cmd'"
					#eval $cmd
#				fi
			done

		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
