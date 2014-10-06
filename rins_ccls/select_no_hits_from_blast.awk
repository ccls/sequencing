BEGIN{
	seq_name=""
	seq_line=""
}
{
	if( $1 == "Query=" ){
		seq_name=$2
		seq_line=$0
	} else if(/No hits found/){
		print seq_name
	}
}
