#!/bin/sh

if [ $# -eq 0 ]; then
	echo "I need at least one filename"
	exit
fi


#	This could be a nice utility script.
#	Import it to sequencing/


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .gz extensino
#	base=${base%.*}	#	drop the unzipped extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \
	srun --share --nice \
		--exclude=n[0000-0019] \
		--job-name="gunzip_${name}" \
		--output=$base \
		--error=$base.gunzip.errors.`date "+%Y%m%d%H%M%S"`  \
		gunzip -c $1 &

#	Sadly, cluster has an old version of gunzip so need to use -c instead of -k

#	gunzip not threaded.  would only block other jobs' usage of node.
#		--cpus-per-task=8 \

#		--output=$base.gunzip.output.`date "+%Y%m%d%H%M%S"`  \
#		gunzip -c $1 > $base &

#	gunzip -c ends up in the srun --output and not my > $base

#	do I need $PWD anywhere? I don't appear to.

	shift
done
