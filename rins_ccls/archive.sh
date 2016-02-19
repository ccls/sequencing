#!/usr/bin/env bash

while [ $# -ne 0 ] ; do
	if [ -f "$1" ] ; then
		echo "Archiving :$1:"
		chmod -f +w md5sums
		chmod -w $1
		md5sum $1 >> md5sums
		gzip --best --verbose $1
		md5sum ${1}.gz >> md5sums
		chmod -w md5sums
	else
		echo "$1 not a file? Skipping"
	fi
shift
done
