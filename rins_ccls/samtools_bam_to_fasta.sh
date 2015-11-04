#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <samfile(s) or bamfile(s)>"
	echo
	echo "Example:"
	echo "  `basename $0` /my/path/*bam"
	echo
	echo "Output is currently a fasta files."
	echo
	echo "Extraction will reverse reads that aligned reverse complemented."
	echo "Extraction will not include reads more that once."
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


while [ $# -ne 0 ] ; do
#	echo $1
	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
	name=${base##*/}	#	just in case given path, drop the path

	#	requires bash >= 4.0
	#	${VARIABLE^^} converts to uppercase
	#	${VARIABLE,,} converts to lowercase
	#if [[ ${ext,,} =~ sam ]] ; then
	#	flag='-S'
	#	samheader="-h $1"
	#else
	#	echo "Input file is not a sam file so not including header in output."
	#	flag=''
	#	samheader=''
	#fi
	[[ ${ext,,} =~ sam ]] && flag='-S -F 256' || flag='-F 256'

#         -f INT   only include reads with all bits set in INT set in FLAG [0]
#         -F INT   only include reads with none of the bits set in INT
#                  set in FLAG [0]
#Flags:
#	0x1	1 PAIRED        .. paired-end (or multiple-segment) sequencing technology
#	0x2	2 PROPER_PAIR   .. each segment properly aligned according to the aligner
#	0x4	4 UNMAP         .. segment unmapped
#	0x8	8 MUNMAP        .. next segment in the template unmapped
#	0x10	16 REVERSE       .. SEQ is reverse complemented
#	0x20	32 MREVERSE      .. SEQ of the next segment in the template is reversed
#	0x40	64 READ1         .. the first segment in the template
#	0x80	128 READ2         .. the last segment in the template
#	0x100	256 SECONDARY     .. secondary alignment
#	0x200	512 QCFAIL        .. not passing quality controls
#	0x400	1024 DUP           .. PCR or optical duplicate
#	0x800	2048 SUPPLEMENTARY .. supplementary alignment

	#	read file
	#	MUST use gawk for the bitwise command "and"
	#	xor() is the same as !and()
	samtools view $flag $1 | gawk '
		(( and( $2 , 4 ) ) || ( xor( $2 , 4 ) && xor( $2, 16 )  )){
			print ">"$1
			print $10
		}
		( xor( $2 , 4 ) &&  and( $2, 16 )  ){
			print ">"$1
			for(i=length($10);i>0;i--) printf substr($10,i,1);
			print " X"
		}' > $base.fasta

#echo welcome | awk '{ for(i=length;i!=0;i--)x=x substr($0,i,1);}END{print x}'
#echo welcome | awk '{ split($0,a,""); for(i=length(a);i>0;printf a[i--]){}}'

	shift
done
