BEGIN{
	if( !num_alignments ) num_alignments=2
	if( !num_descriptions ) num_descriptions=2
}
{
	if( $1 == "Query=" ){
		align_count=0
		desc_count=-1
	}
	if( /^>/ ){
		align_count++
		desc_count=-1
	}
	if( /^Lambda     K      H/ ){
		align_count=0
	}


	if( desc_count >= 0 ){
		if( !/^\s*$/ ){
			desc_count++
		}
	}
	if( desc_count > 3 ){
		if( /^\s*$/ ){
			desc_count=-1
		}
	}
	if( /^Sequences producing significant alignments/ ){
		desc_count=0
	}


	if( ( align_count <= num_alignments ) && ( desc_count <= num_descriptions ) ){
		print $0
	}
}
