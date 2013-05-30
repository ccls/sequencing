#!/bin/sh

if [ `uname -n` = "fxdgroup-169-229-196-225.sph.berkeley.edu" ] ; then
	#	echo "on dev"
	db='/Volumes/cube/working/indexes/nt'
	cmdbase=''
else
	db='/my/home/jwendt/blast/nt'
	cmdbase='srun'	#	pointless when only running one ... --exclusive -n 1'
fi

#	each filename on the command line
while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then
		now=`date "+%Y%m%d%H%M%S"`
		subdir=$1.${now}.pieces
		mkdir $subdir

#	different version of split on ec !!!!
#		split -a 4 -p '>' $1 $subdir/${1}_

		#	This seems to do the trick
		awk 'BEGIN{c=0}{if(/^>/){f=sprintf("'$subdir/$1'_%09d",++c)}print>>f}' $1

#
#	some code suggests that we may need to close a file in awk????
#

		for file in `ls $subdir/${1}_*` ; do
			cmd="$cmdbase blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt &"
			echo $cmd
			#	need to eval to use the &
		#	eval $cmd
		done

	else
		echo "$1 doesn't seem to be a file."
	fi
	shift
done #	while [ $# -ne 0 ] ; do
