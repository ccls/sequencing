BEGIN{
	if( head == "" ) head="BLASTN"
	if( tail == "" ) tail="Gap Penalties"
	head_count=0
	tail_count=0
}
(NR%100000 == 0){
	print "Read",FNR,"lines from",FILENAME >> "blast_check.log"
}
( $0 ~ head ){ head_count++ }
( $0 ~ tail ){ tail_count++ }
(/Query= /){ query++ }
(/Effective search space used:/){ effective++ }
(/no longer exists in database/){ nolonger++ }
(/[^[:print:]]/){ 	#	too many (includes \anything?)
	print "Non-printable character(s) found on line number :",NR,":"
	print $0
	nonprint++ 
}
(/[[:cntrl:]]/){	#	REQUIRES double brackets as they are for different reasons
	print "Control character(s) found on line number :",NR,":"
	print $0
	control++ 
}
END{
	print "First line count :",head_count,":"
	print "Last line count :",tail_count,":"
	if( head_count != tail_count ){ print " * First and Lines lines are out of sync" }
	print "Query=  line count :",query,":"
	print "Effective search space used: line count :",effective,":"
	if( query != effective ){ print " * Query= and Effective search lines are out of sync" }
	print "no longer exists in database line count :",nolonger,":"
	if( nolonger > 0 ){ print " * TOO MANY no longer exists in database lines" }
	print "nonprint character line count :",nonprint,":"
	if( nonprint > 0 ){ print " * TOO MANY nonprintable characters" }
	print "control character line count :",control,":"
	if( control > 0 ){ print " * TOO MANY control characters" }
	print "---"
	close("blast_check.log")
}
