#!/usr/bin/env bash

#	50
#	Start of HERV K113 LTR
#ltr_start="TGTGGGGAAAAGCAAGAGAGATCAGATTGTTACTGTGTCTGTGTAGAAAG"
#	End of HERV K113 LTR
#ltr_end="CCACCTTACGAGAAACACCCACAGGTGTGTAGGGGCAACCCACCCCTACA"

#	75
#	Start of HERV K113 LTR
#ltr_start="TGTGGGGAAAAGCAAGAGAGATCAGATTGTTACTGTGTCTGTGTAGAAAGAAGTAGACATAGGAGACTCCATTTT"
#	End of HERV K113 LTR
#ltr_end="CTTTTTCTTTTCCAAATCTCTCGTCCCACCTTACGAGAAACACCCACAGGTGTGTAGGGGCAACCCACCCCTACA"

#	90
#	Start of HERV K113 LTR
ltr_start="TGTGGGGAAAAGCAAGAGAGATCAGATTGTTACTGTGTCTGTGTAGAAAGAAGTAGACATAGGAGACTCCATTTTGTTATGTACTAAGAA"
#	End of HERV K113 LTR
ltr_end="ACTTTGTCTCTGTGTCTTTTTCTTTTCCAAATCTCTCGTCCCACCTTACGAGAAACACCCACAGGTGTGTAGGGGCAACCCACCCCTACA"

#for line in `awk -F, '{print $1}' overlappers_table.Q20.csv | tail +2` ; do
#	different versions of tail.  gnu tail requires a -n
#for line in `awk -F, '{print $1}' overlappers_table.ALL.csv | tail -n +2` ; do
#for line in `awk -F, '{print $1}' grouping.overlappers_table.ALL.csv | tail -n +2` ; do
for line in `awk -F, '{print $1}' overlappers_table.Q20.csv | tail -n +2` ; do

#       ${parameter#word}
#       ${parameter##word}
#              Remove matching prefix pattern.
#       ${parameter%word}
#       ${parameter%%word}
#              Remove  matching  suffix  pattern.  

#	echo $line
#chrY:8311408:R
#	echo ${line#*:}
#8311408:R
#	echo ${line##*:}
#R
	direction=${line##*:}
#	echo ${line%:*}
#chrY:8311408
	insertion_point=${line%:*}
	position=${insertion_point##*:}
#	echo ${line%%:*}
#chrY
	chromosome=${line%%:*}

#	echo $chromosome
#	echo $position
#	echo $insertion_point
#	echo $direction	#	F or R


#	CONSIDER DIRECTION!!!!!
#	Does it really matter in this context?


#	if [[ $direction == 'F' ]] ; then

#	using the ,, or ^^ requires bash > 4.0
# requires bash >= 4.0
# ${VARIABLE^^} converts to uppercase
# ${VARIABLE,,} converts to lowercase
###if [[ ${ext,,} =~ sam ]] ; then


		#	The insertion point is the first position of the ltr so DO NOT INCLUDE
		let to=position-1
		let from=to-89	#74	#49
		echo ">${chromosome}:${from}_${to}:${direction}+HERVK113STARTLTR"
		hg=`samtools faidx FASTA/hg19.fa ${chromosome}:${from}-${to} | grep -vs "^>" | tr -d '\n'`
		echo ${hg^^}
		echo $ltr_start

		#	expecting overlap of 6 (although not always 6)
		
		let from=position-7
		let to=from+89	#74	#49
		echo ">HERVK113ENDLTR+${chromosome}:${from}_${to}:${direction}"
		echo $ltr_end
		hg=`samtools faidx FASTA/hg19.fa ${chromosome}:${from}-${to} | grep -vs "^>" | tr -d '\n'`
		echo ${hg^^}

#	else

#	I must understand what it really means when a read matches as reverse complement
#	I see a read that says it matches RC, but when I look at the raw read and compare
#	it to whats in the bam file and the genome, they look normal.  Confused.
#
#	@ERR251139.54934541 FCD1LREACXX:6:2316:16546:77353/2
#	TTCACCCTAGAGAAAAGCCTCCACGTTGGGCACCAGA TGTAGGGGTGGGTTGCCCCTACACACCTGTGGGTGTTTCTCGTAAGGTGGGACGAGAGGTTTG
# |- Normal Match yet flagged RC???-->| |<--- Reverse Complement of HERV END LTR    <-------
#                                       TGTAGGGGTGGGTTGCCCCTACACACCTGTGGGTGTTTCTCGTAAGGTGG
#
#	There is a command "rev" that reverses the order of a string
#	that could be useful. "echo 1234 | rev" -> "4321"

#	Create the reverse complement like so ...
#echo "TGTGGGGAAAAGCAAGAGAGATCAGATTGTTACTGTGTCTGTGTAGAAAG" | rev | tr "ATCG" "TAGC"
#CTTTCTACACAGACACAGTAACAATCTGATCTCTCTTGCTTTTCCCCACA

#	fi


done
