#!/usr/bin/env bash

if [ $# -lt 1 ]; then
	echo "I need at least 1 fast[qa] filename."
	exit
fi

base=`basename $PWD`

srun --nice --share --partition=all \
	--exclude=n[0000-0009] \
	--job-name="overlapper_seeker_${base}" \
	--cpus-per-task=8 \
	--error=$base.overlapper_seeker.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
	--output=$base.overlapper_seeker.output.`date "+%Y%m%d%H%M%S"`.nobackup \
	overlapper_seeker.sh $* &

