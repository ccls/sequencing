#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--other OTHER_BAM2FASTX_OPTIONS] <bamfile(s)>"
	echo
	echo "Example:"
	echo "  `basename $0` --other '-Q --paired' /my/path/*bam"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

other=''

while [ $# -ne 0 ] ; do
	case $1 in
		-o|--o*)
			shift; other=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .bam extension
	name=${base#*/}	#	just in case given path

	cmd="bam2fastx $other --fastq --all -N -o $name.fastq $1"
	echo $cmd

#		--begin=23:00 \
#		--partition=bigmem \
#		--exclude=n[0000-0029] \

	srun --nice --share \
		--partition=all \
		--exclude=n[0000-0009] \
		--job-name="bam2fastx_${name}" \
		--cpus-per-task=8 \
		--error=$base.bam2fastx.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.bam2fastx.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		$cmd &

#	bam2fastx is from tophat

	shift
done
