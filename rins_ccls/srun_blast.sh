#!/bin/sh

function usage(){
	echo
	echo "Creates a directory of the name of each of the given fasta files and appends"
	echo "the date and .pieces.  The fasta file is then split into fasta files each"
	echo "with no more than the maximum number of reads and sequentially numbers them"
	echo "placing them in the aforementioned directories."
	echo
	echo "Usage: (NO EQUALS SIGNS)"
	echo
	echo "`basename $0` [--command COMMAND] [--db BLAST_DATABASE] "
	echo "              [--evalue FLOAT]"
	echo "              [--outfmt BLASTN_OUTFMT#] fasta_filelist"
	echo
	echo "Defaults:"
	echo "  command . : blastn"
	echo "  db  ..... : nt"
	echo "  outfmt .. : 0"
	echo "  ......... :   0=pairwise text"
	echo "  ......... :  10=csv"
	echo "  evalue .. : 0.05"
	echo 
	echo "Example:"
	echo "  `basename $0` --db nt /my/path/*fasta"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

outfmt=0
evalue=0.05
command='blastn'
db='nt'

while [ $# -ne 0 ] ; do
	case $1 in
		-c|--c*)
			shift; command=$1; shift ;;
		-d|--d*)
			shift; db=$1; shift ;;
		-e|--e*)
			shift; evalue=$1; shift ;;
		-o|--o*)
			shift; outfmt=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done


while [ $# -ne 0 ] ; do
	echo $1
	base=${1%.*}		#	drop the .fasta extension
	name=${base#*/}	#	just in case given path

#		--begin=23:00 \
#		--partition=bigmem \

	srun --share --nice \
		--partition=all \
		--exclude=n[0000-0009] \
		--begin=23:00 \
		--job-name="${command}_${db}_${name}" \
		--cpus-per-task=8 \
		--output=$base.${command}_${db}.${evalue}.output.`date "+%Y%m%d%H%M%S"`  \
		--error=$base.${command}_${db}.${evalue}.errors.`date "+%Y%m%d%H%M%S"`  \
		${command} -num_threads 8 -num_alignments 20 -num_descriptions 20 \
			-evalue ${evalue} -outfmt ${outfmt} -db ${db} \
			-query $1 \
			-out $base.${command}_${db}.${evalue}.txt &

	shift
done
