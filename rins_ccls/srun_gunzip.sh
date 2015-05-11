#!/usr/bin/env bash

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .gz extensino
#	base=${base%.*}	#	drop the unzipped extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \
	srun --share --nice \
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="gunzip_${name}" \
		--output=$base.gunzip.output.`date "+%Y%m%d%H%M%S"`.nobackup  \
		--error=$base.gunzip.errors.`date "+%Y%m%d%H%M%S"`.nobackup  \
		gunzip -k $1 &

#	20150317 - gzip upgraded! Can use -k now. --output should be empty now.
#		--output=$base \

#	Sadly, cluster has an old version of gunzip so need to use -c instead of -k

#	gunzip not threaded.  would only block other jobs' usage of node.
#		--cpus-per-task=8 \

#	gunzip -c ends up in the srun --output and not my > $base

	shift
done
