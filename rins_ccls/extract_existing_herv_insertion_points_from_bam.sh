#!/usr/bin/env bash

#
#	Given a bam file containing the alignments
#	to HERV beginning and ending LTR
#	this script will extract the insertion point.
#	(the beginning of the beginning or the ending or the ending LTR)
#

#	B = beginning
#	E = ending
#	F = forward
#	R = reverse

while [ $# -ne 0 ] ; do

	samtools view -F 20 $1 | grep "beginning" | awk '{print $3":"$4":BF"}' | sort
	samtools view -F 20 $1 | grep "ending" | awk '{print $3":"$4+length($10)":EF"}' | sort

	samtools view -F 4 -f 16 $1 | grep "beginning" | awk '{print $3":"$4+length($10)":BR"}' | sort
	samtools view -F 4 -f 16 $1 | grep "ending" | awk '{print $3":"$4":ER"}' | sort

	shift
done
