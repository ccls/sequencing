#!/bin/sh

#if [ `uname -n` = "fxdgroup-169-229-196-225.sph.berkeley.edu" ] ; then
	#	echo "on dev"
	db='/Volumes/cube/working/indexes/nt'
	cmdbase=''
#else
#if [ `uname -n` = "genepi1.berkeley.edu" ] ; then
if [ `uname -n` = "ec0000" ] ; then
	db='/my/home/jwendt/dna/blast/nt'
	cmdbase='simple_queue.sh push srun --share'	#	pointless when only running one ... --exclusive -n 1'
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
		fasta_file=`basename $1`
#	not used so why bother
#		fasta_dir=`dirname $1`
		subdir=$fasta_file.${now}.pieces
		mkdir $subdir

		awk '
			BEGIN{
				file_number=0
				read_count=0
				f=sprintf("'$subdir/$fasta_file'_%09d",++file_number)
			}
			{
				if(/^>/){
					if( read_count >= '$max_reads' ){
						close(f)
						f=sprintf("'$subdir/$fasta_file'_%09d",++file_number)
						read_count=0
					}
					read_count++
				}
				print>>f
			}' $1

		for file in `ls $subdir/${fasta_file}_*` ; do
			cmd="$cmdbase blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt &"
			echo $cmd

#			if [ `uname -n` = "genepi1.berkeley.edu" ] ; then
			if [ `uname -n` = "ec0000" ] ; then
				#	need to eval to use the &
				eval $cmd
			fi

		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
