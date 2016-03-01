#
#	Assuming that file contains ...
#		position,TCGA-02-2483-01A,TCGA-02-2483-10A,...
#		chr1:1345187:BF,29,25,..
#		chr1:1346154:EF,17,16,..
#	... and is sorted like "sort -t: -k1,1 -k2,2n"
#
#	Note: all matching blocks run for each line
#
BEGIN {
	FS=","
	OFS=","
}
function set_buffer() {
	bufferer_position[1]=position[1]
	bufferer_position[2]=position[2]
	for(i=1;i<=NF;i++){b[i]=$i}
}
function print_buffer(){
	out=""
	for(i=1;i<NF;i++){out=out""b[i]","}
	print out""b[NF]
}
( NR == 1 ){ print }
( NR > 1 ){
	split($1,position,":")
	if( bufferer_position[1] ){
		if(( position[1] == bufferer_position[1] ) && ( position[2] < (bufferer_position[2]+10000) )){
			for(i=2;i<=NF;i++){b[i]+=$i}
		} else {
			print_buffer();
			set_buffer();
		}
	} else {
		set_buffer();
	}
}
END {
	print_buffer();
}
