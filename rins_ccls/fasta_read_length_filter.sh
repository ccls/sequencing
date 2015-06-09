#!/usr/bin/env bash

function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` [--length INTEGER] [--side l/r] <fasta files(s)>"
	echo
	echo "Example:"
	echo "  `basename $0` -s r -l 30 /my/path/*fasta"
	echo
	echo "Defaults:"
	echo "  --length/-l ..... : 50"
	echo "  --side/-s   ..... : l"
	echo
	echo "Output directed to file, similarly named."
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage

length=50
side='l'

while [ $# -ne 0 ] ; do
	case $1 in
		-l|--l*)
			shift; length=$1; shift ;;
		-s|--s*)
			shift; side=$1; shift ;;
		-*)
			echo ; echo "Unexpected args from: ${*}"; usage ;;
		*) 
			break;;
	esac
done

while [ $# -ne 0 ] ; do
#	echo $1
	base=${1%.*}		#	drop the extension
	ext=${1##*.}		#	grab the extension
#	name=${base#*/}	#	just in case given path, drop the path

	awk -v l=$length -v s=$side -v out=$base.trim$length.$ext 'BEGIN{ n=r="" }
		END{ o() }
		function o(){
			if( length(r) >= l ){
				print n >> out
				if( s == "l") {
					print substr(r,1,l) >> out
				} else {
					print substr(r,length(r)-l+1,l) >> out
				}
			}
		}
		( /^>/ ){
			o()
			n=$0
			r=""
		}
		( !/^>/ ){ r=r $0 }' $1

	shift
done
