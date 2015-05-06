#!/bin/sh

function usage(){
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

outfmt='0'
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

#[[ $outfmt -eq '10' ]] && ext='csv' || ext='txt'
#	-eq is for numbers.  == is for strings.
[[ $outfmt == '10' ]] && ext='csv' || ext='txt'

#	Warning: The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output

if [[ $outfmt -gt 4 ]] ; then
	other='-max_target_seqs 20'
else
	other='-num_alignments 20 -num_descriptions 20'
fi

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
		--output=$base.${command}_${db}.${evalue}.output.`date "+%Y%m%d%H%M%S"`.nobackup  \
		--error=$base.${command}_${db}.${evalue}.errors.`date "+%Y%m%d%H%M%S"`.nobackup  \
		${command} -num_threads 8 $other \
			-evalue ${evalue} -outfmt ${outfmt} -db ${db} \
			-query $1 \
			-out $base.${command}_${db}.${evalue}.${ext} &

	shift
done
