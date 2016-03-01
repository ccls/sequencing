#	if R2 and seq name is same as buffered, print both
#	MUST use gawk for the bitwise command "and"
#	samtools view $flag $1 | gawk '

#	THIS WAS JUST SNIPPED FROM A SCRIPT.
#	IT NEED ${base} FIXED 

( and( $2 , 64 ) ){
	b1=$1
	b10=$10
	b11=$11
}
( and( $2 , 128 ) ){
	if ( $1 == b1 ){
		if ( length(b10) != length($10) ){
			print $1 >> "'${base}.diff_length_reads'"
		}
		if ( length(b11) != length($11) ){
			print $1 >> "'${base}.diff_length_quality'"
		}
		print "@"b1"/1" >> "'${base}_R1.fastq'"
		print b10 >> "'${base}_R1.fastq'"
		print "+" >> "'${base}_R1.fastq'"
		print b11 >> "'${base}_R1.fastq'"

		print "@"$1"/2" >> "'${base}_R2.fastq'"
		print $10 >> "'${base}_R2.fastq'"
		print "+" >> "'${base}_R2.fastq'"
		print $11 >> "'${base}_R2.fastq'"
		b1=b10=b11=""
	}
}
