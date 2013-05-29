#!/bin/sh

if [ `uname -n` = "fxdgroup-169-229-196-225.sph.berkeley.edu" ] ; then
	#	echo "on dev"
	db='/Volumes/cube/working/indexes/nt'
	cmdbase=''
else
	db='/my/home/jwendt/blast/nt'
	cmdbase='srun --exclusive -n 1'
fi

while [ $# -ne 0 ] ; do

	if [ -f $1 ] ; then
		cmd="$cmdbase blastn -query $1 -db $db -evalue 0.05 -outfmt 0 -out $1.blastn.txt &"
		touch $1.blastn.txt.start_touch
		eval $cmd	#	need to eval to use the &
	else
		echo "$1 doesn't seem to be a file."
	fi

	shift
done #	while [ $# -ne 0 ] ; do

