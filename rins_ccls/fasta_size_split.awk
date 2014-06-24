#!/usr/bin/env awk -f
#
##!/bin/sh
##
##	split fasta files into smaller fasta files with a 100000+ characters
##
##		fasta_size_split.sh my_file_1.fasta my_file_2.fasta
##
##	if first argument is a number, it will be used for the size
##
#
#
#if [ $# -eq 0 ]; then
#	echo
#	echo "split fasta files into smaller fasta files with a 100000+ characters"
#	echo
#	echo "Usage:"
#	echo
#	echo "`basename $0` [optional max_characters_ish]"
#	echo
#	echo "Example:"
#	echo
#	exit
#fi
#
#tmp=`echo $1 | tr -cd '[:digit:]'`
## need the x's in case is blank
#if [ "x${tmp}" == "x${1}" ] ; then
#	size=$1
#	shift
#else
#	size=100000
#fi
#
#while [ $# -ne 0 ] ; do
#	awk '
	BEGIN{
		size=100000	#	explicitly overridable on the command line
		current_file_name=""
	}
	( FILENAME != current_file_name ){ 
		current_file_name=FILENAME
		file_number=0
		char_count=0
		f=sprintf("%s_%04d",FILENAME,++file_number)
	}
	{
		char_count+=length
		if(/^>/ && char_count >= size){
			close(f)
			f=sprintf("%s_%04d",FILENAME,++file_number)
			char_count=0
		}
		print>>f
	}
#	' $1
#		shift
#	done
