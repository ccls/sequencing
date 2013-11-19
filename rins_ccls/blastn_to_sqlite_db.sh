#!/bin/sh

if [ $# -eq 0 ]; then
	echo
	echo "No files given"
	echo
	echo "Usage:"
	echo
	echo "$0"
	echo
	echo "Example:"
	echo "$0 trinity_input_paired.blastn.txt"
	echo
	exit
fi

#
#Query= comp14_c0_seq1 len=108 path=[1:0-46 48:47-107]
#
#Length=108
#                                                                      Score     E
#Sequences producing significant alignments:                          (Bits)  Value
#
#gb|AC027811.9|  Homo sapiens chromosome 15, clone RP11-358L4, com...   158    1e-35
#ref|NG_011635.1|  Homo sapiens myosin IIIA (MYO3A), RefSeqGene on...   119    4e-24
#emb|AL391812.7|  Human DNA sequence from clone RP13-94I22 on chro...   119    4e-24
#emb|AL590407.8|  Human DNA sequence from clone RP11-522L3 on chro...   115    6e-23
#gb|AC114279.2|  Homo sapiens chromosome 5 clone CTD-2281M20, comp...   113    2e-22
#emb|AL365189.15|  Human DNA sequence from clone RP11-595C20 on ch...   106    3e-20
#gb|AC147049.2|  Pan troglodytes BAC clone RP43-131M1 from chromos...  99.0    6e-18
#gb|AC073275.8|  Homo sapiens BAC clone RP11-651K8 from 7, complet...  99.0    6e-18
#gb|AC200162.4|  Pan troglodytes BAC clone CH251-643F16 from chrom...  93.5    3e-16
#gb|AC146195.3|  Pan troglodytes BAC clone RP43-25F5 from chromoso...  93.5    3e-16
#gb|AC099798.4|  Homo sapiens BAC clone RP11-559E21 from 7, comple...  93.5    3e-16
#gb|AC073181.8|AC073181  Homo sapiens clone RP13-493F23, complete ...  93.5    3e-16
#emb|AL357872.23|  Human DNA sequence from clone RP11-308P9 on chr...  86.1    5e-14
#
#
#>gb|AC027811.9| Homo sapiens chromosome 15, clone RP11-358L4, complete sequence
#Length=175664
#
# Score =  158 bits (85),  Expect = 1e-35
# Identities = 87/88 (99%), Gaps = 0/88 (0%)
# Strand=Plus/Plus
#
#Query  1      GAAAGAAAGATACAGTCTTTTTCAGATGAACAAATGCTGAGAGAATTCGCCACTCCCAGG  60
#              ||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#Sbjct  35260  GAAGGAAAGATACAGTCTTTTTCAGATGAACAAATGCTGAGAGAATTCGCCACTCCCAGG  35319
#
#Query  61     CTACCACTACAAGAACTGCCAAAAGGAG  88
#              ||||||||||||||||||||||||||||
#Sbjct  35320  CTACCACTACAAGAACTGCCAAAAGGAG  35347
#
#
#

while [ $# -ne 0 ] ; do
	if [ -f $1 ] ; then

		datestamp=`date "+%Y%m%d%H%M%S"`
		db="$1.$datestamp.sqlite3"

		sqlite3 -cmd '.timeout 5000' $db 'create table blastn(id integer primary key autoincrement, query, len, result, bitscore, score, expect, identities, identities_percent, gaps, gaps_percent, strand)'

		#	can't use "length" in awk. built-in function

		#Length= is used twice and in different ways
		#	want the length of the query, NOT the length of the entire aligned sequence
		#	so only set it if its blank so we catch the first one

		#	result (the matched sequence name) is only the first line in the blastn output
		#		Sometimes, this name is multiple lines.  Is it worth the trouble to get the whole thing?

		#	using octal codes for quotes
		#	\42 = "
		#	\47 = '

		#	removing apostrophes in dq function
		#	can't seem to find a way to escape them
		#	Children's Hospital
		#	'\''

		#	this works
		#			gsub("\47","\47\134\47\47",value)
		#			gsub("\x27","\x27\x5C\x27\x27",value)
		#	can't use single quotes in the awk program, but can use the octal (or hex) codes
		#	escape the middle quote and voila

		#	double up a double quote

		awk '
		BEGIN { reset() }
		function reset(){
			query=""
			len=""
			result=""
			bitscore=""
			score=""
			expect=""
			identities=""
			identities_percent=""
			gaps=""
			gaps_percent=""
			strand=""
		}
		function dq(value){
			gsub("\47","\47\134\47\47",value)
			gsub("\42","\42\42",value)
			return "\42"value"\42"
		}
		function show(){
			print result
			print len, bitscore, score, expect, identities, identities_percent, gaps, gaps_percent, strand

			command="sqlite3 -cmd \47.timeout 5000\47 '$db' \47insert into blastn(" \
				"query,result,len,bitscore,score,expect,identities,identities_percent," \
				"gaps,gaps_percent,strand) values (" dq(query) "," dq(result) "," \
					len "," \
					bitscore "," \
					score "," \
					dq(expect) "," \
					dq(identities) "," \
					identities_percent "," \
					dq(gaps) "," \
					gaps_percent "," \
					dq(strand) ")\47"
			print command
			status=system(command)
			if( status != 0 ){ exit }
		}
		{
			if( $1 == "Query=" ){
				if( result != "" ){ show() }
				reset()
				query=$2
				print "\n" $0
			} else if( /Length=\d*/ ){
				if( len == "" ){ len=substr($1,8) }
			} else if( /^>/ ){
				if( result != "" ){ show() }
				result=substr($0,2)
			} else if( $1 == "Score" ){
				bitscore=$3
				score=$5
				gsub(/\(|\)|,/,"",score)
				expect=$8
			} else if( $1 == "Identities" ){
				identities=$3
				identities_percent=$4
				gsub(/\(|%|,|\)/,"",identities_percent)
				gaps=$7
				gaps_percent=$8
				gsub(/\(|%|\)/,"",gaps_percent)
			} else if( /Strand=/ ){
				strand=substr($1,8)
			}
		}' $1
	else
		echo "File '$1' not found?"
	fi
	shift
done

