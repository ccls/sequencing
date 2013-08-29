#!/bin/sh


uname=`uname -n`

db='/Volumes/cube/working/indexes/nt'

if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
	db='/my/home/jwendt/dna/blast/nt'
fi

tmp=`echo $1 | tr -cd '[:digit:]'`
#	need the x's in case is blank
if [ "x${tmp}" = "x${1}" ] ; then
	max_reads=$1
	shift
else
	max_reads=1000
fi

#	each filename on the command line
while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then
		now=`date "+%Y%m%d%H%M%S"`

#		fasta_base=`basename $1`
		fasta_base=`basename $1 | sed 's/\(.*\)\..*/\1/'`

		#	extension=${filename##*.}
		#	the longest matching  pattern  (the  ``##'' case)  deleted
		#	It deletes the longest *. match leaving only the extension behind (starts at the beginning (left))
		#	% and %% start from the trailing or right side (could be similar to using basename)

		#	not used so why bother
		#		fasta_dir=`dirname $1`

#		subdir=$fasta_base.${now}.pieces
		subdir=`basename $1`.${now}.pieces

		#
		#	if given path to file, put output there or where command was executed???
		#	I say where command is executed.
		#
		mkdir $subdir

		#
		#	I would like to drop the original fasta file's extension,
		#	but that would require changes below.  Be careful.
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
		for file in `find $PWD/$subdir/ -type f -name ${fasta_base}.\*.fasta` ; do
			cmd=''	#	gotta reset it
			if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
#	trinity_input_single.uniq.fasta_00000416.fasta
#	=> 'uniq' for num.  oops
#				num=`basename $file | awk -F. '{print $2}' | awk -F_ '{print $NF}'`
#				num=`basename $file | awk -F. '{print $(NF-1)}' | awk -F_ '{print $NF}'`
				num=`basename $file | awk -F. '{print $(NF-1)}'`
				cmd="srun --share --job-name=$num"
			fi
			cmd="$cmd blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt &"


#			cmd="$cmd blastn_wrapper.sh $file $db &"

			echo $cmd

			if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
				#	need to eval to use the &
				#	want the & in the queue'd command, not here.
				#	if were here, will cause database error by trying to write to it at same time
				eval "simple_queue.sh push '$cmd'"
				#eval $cmd
			fi

		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
