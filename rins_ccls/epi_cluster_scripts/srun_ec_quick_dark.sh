#!/usr/bin/env bash

if [ $# -ne 2 ]; then
	echo "I need EXACTLY 2 filenames.  1 left lane (R1).  1 right lane (R2)."
	exit
fi

base=`basename $PWD`

srun --nice --share --partition=bigmem \
	--job-name="dark_${base}" \
	--cpus-per-task=8 \
	--error=$base.quick_dark.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
	--output=$base.quick_dark.output.`date "+%Y%m%d%H%M%S"`.nobackup \
	ec_quick_dark.sh $1 $2 &
