#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo
	echo "No files given"
	echo
	echo "Usage:"
	echo
	echo "Groups and sums table rows if within 10k bp of previous."
	echo "(this can get a little wierd if row 2 is 9k from row 1 and row 3 is 9k from row 2)"
	echo "(row 1 and 2 would get grouped and row 3 would be on its own)"
	echo "(grouping won't go backwards)"
	echo
	echo "`basename $0` csv_table"
	echo
	echo "position,TCGA-02-2483-01A,TCGA-02-2483-10A,TCGA-02-2485-01A,TCGA-02-2485-10A"
	echo "chr1:1345187:BF,76,58,65,51"
	echo "chr1:10487684:EF,1,1,,2"
	echo "chr1:11365577:BR,5,1,2,2"
	echo "chr1:12068255:BR,4,4,5,2"
	echo "chr1:15462793:BR,19,7,14,14"
	echo
	echo "Example:"
	echo "$0 insertion_points_near_reference.hg19.Q20.csv"
	echo
	exit
fi

while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then

		dir=`dirname $0`
		awk -f "$dir/csv_table_group_rows.awk" $1

	else
		echo "File '$1' not found?"
	fi
	shift
done
