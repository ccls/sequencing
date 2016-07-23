#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo
	echo "No files given"
	echo
	echo "Usage:"
	echo
	echo "`basename $0` 3_column_csv_file"
	echo
	echo "chr13:20175326:F:POST,2,HG00096"
	echo "chr14:38588273:F:POST,2,HG00096"
	echo "chr15:28430105:F:POST,1,HG00096"
	echo "chr15:66026960:F:POST,1,HG00096"
	echo "chr15:89084751:F:POST,2,HG00096"
	echo "chr16:34234143:F:POST,2,HG00096"
	echo
	echo "Example:"
	echo "$0 insertion_points.hg19.Q20"
	echo
	exit
fi

while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then

		dir=`dirname $0`
		#	gawk is needed as this program uses multidimensional arrays
		gawk -f "$dir/to_table.gawk" $1

	else
		echo "File '$1' not found?"
	fi
	shift
done
