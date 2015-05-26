#!/usr/bin/env bash

if [ $# -lt 1 ]; then
	echo "I need at least 1 fast[qa] filename."
	exit
fi

base=`basename $PWD`

srun --nice --share --partition=all \
	--exclude=n[0000-0009] \
	--job-name="herv_k113_${base}" \
	--cpus-per-task=8 \
	--error=$base.herv_k113.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
	--output=$base.herv_k113.output.`date "+%Y%m%d%H%M%S"`.nobackup \
	herv_k113_integration.sh $* &
