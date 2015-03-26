#!/bin/sh

function usage(){
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--db BOWTIE2_INDEX] "
	echo
	echo "Defaults:"
	echo "  db  ..... : hg19"
	echo 
	echo "Example:"
	echo "  `basename $0` --db hg19 /my/path/*fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

db='hg19'

while [ $# -ne 0 ] ; do
	case $1 in
		-d|--d*)
			shift; db=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done

while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the extension
	name=${base#*/}	#	just in case given path

	#	This better be a or q
	#filetype=${1:(-1)}
	[[ ${1:(-1)} -eq 'q' ]] && filetype='-q' || filetype='-f'


#		--partition=bigmem \
#		--exclude=n[0000-0029] \
#		--begin=23:00 \

	srun --nice --share --partition=bigmem \
		--job-name="bowtie2_${name}_${db}" \
		--cpus-per-task=8 \
		--error=$base.bowtie2.${db}.very-sensitive-local.errors.`date "+%Y%m%d%H%M%S"`.nobackup \
		--output=$base.bowtie2.${db}.very-sensitive-local.output.`date "+%Y%m%d%H%M%S"`.nobackup \
		bowtie2 $filetype --threads 8 \
			--very-sensitive-local \
			-x $db \
			-U $1 -S $base.bowtie2.${db}.very-sensitive-local.sam &

#	bowtie2 can use $BOWTIE2_INDEXES for path
#			-x /my/home/ccls/indexes/bowtie2/herv_k113 \

#  For --local:
#   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)

#	output goes to STDERR (--error)

	shift
done
