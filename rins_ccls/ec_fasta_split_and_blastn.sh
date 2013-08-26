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

#
#	if given path to file, put output there or where command was executed???
#	I say where command is executed.
#
		fasta_base=`basename $1`
#	not used so why bother
#		fasta_dir=`dirname $1`
		subdir=$fasta_base.${now}.pieces
		mkdir $subdir

		awk '
			BEGIN{
				file_number=0
				read_count=0
				f=sprintf("'$subdir/$fasta_base'_%08d.fasta",++file_number)
			}
			{
				if(/^>/){
					if( read_count >= '$max_reads' ){
						close(f)
						f=sprintf("'$subdir/$fasta_base'_%08d.fasta",++file_number)
						read_count=0
					}
					read_count++
				}
				print>>f
			}' $1

		#	on some occassions, this list is too long for ls so changing to find
		#	for file in `ls $PWD/$subdir/${fasta_base}_*.fasta` ; do
		for file in `find $PWD/$subdir/ -type f -name ${fasta_base}_*.fasta` ; do
			cmd=''	#	gotta reset it
			if [ $uname = "ec0000" -o $uname = "n0.berkeley.edu" ] ; then
#	trinity_input_single.uniq.fasta_00000416.fasta
#	=> 'uniq' for num.  oops
#				num=`basename $file | awk -F. '{print $2}' | awk -F_ '{print $NF}'`
				num=`basename $file | awk -F. '{print $(NF-1)}' | awk -F_ '{print $NF}'`
				cmd="srun --share --job-name=$num"
			fi
			cmd="$cmd blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt &"

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
