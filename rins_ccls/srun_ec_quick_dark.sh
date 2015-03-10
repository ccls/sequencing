#!/bin/sh

if [ $# -ne 2 ]; then
	echo "I need EXACTLY 2 filenames.  1 left lane (R1).  1 right lane (R2)."
	exit
fi

base=`basename $PWD`

#while [ $# -ne 0 ] ; do
#	echo $1
#	base=${1%.*}		#	drop the .fasta extension
#	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

srun --nice --share --partition=bigmem \
	--job-name="quick_dark${base}" \
	--cpus-per-task=8 \
	--error=$base.quick_dark.errors.`date "+%Y%m%d%H%M%S"` \
	--output=$base.quick_dark.output.`date "+%Y%m%d%H%M%S"` \
	ec_quick_dark.sh $1 $2 &
