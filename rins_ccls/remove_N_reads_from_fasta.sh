#!/bin/sh

#
#	exec or no exec.  Here, really makes no difference.
#
#	sequence may be multiple lines.  $0 trims newline.
#	keeping newline in concatenation, so need to use
#	printf when printing it as print will add one.
#

#	I don't see the purpose of using exec here.
#	May be none as is last statement in script.
#exec awk '

awk '
function maybe_print(){
	if( ( sequence!="" ) && ( sequence !~ /N/ ) ){
		print name
		printf sequence
	}
}
/^>/{ 
	maybe_print()
	name=$0 
	sequence=""
} 
!/^>/{ 
	sequence=sequence$0"\n"
}
END{
	maybe_print()
}' $@


#	ONLY WORKS WHEN SEQUENCE IS ALL ON ONE LINE!!!!
#!/(^>)|N/{
#	print name
#	print $0
#}' $@
