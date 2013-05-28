#!/bin/sh

db='/Volumes/cube/working/indexes/nt'
#	-db=/my/home/jwendt/blast/nt \

#	gonna try to use this on the Epi Cluster

while [ $# -ne 0 ] ; do
	echo $1

	if [ -f $1 ] ; then
		echo "blastn -query=$1 -db=$db -evalue 0.05 -outfmt 0 > $1.blastn.txt"
	else
		echo "$1 doesn't seem to be a file."
	fi

	shift
done #	while [ $# -ne 0 ] ; do
