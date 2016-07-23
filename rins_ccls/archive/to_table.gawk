BEGIN{
	FS=","
}
{
	p[$1]++
	s[$3]++
	b[$1][$3]=$2
}
END{
	asorti(p)
	asorti(s)
	printf "position"
	for(subj in s)
		printf ",%s",s[subj]
	printf "\n"

	for(pos in p){
#		sum=0;
#		for(subj in s)
#			sum+=b[p[pos]][s[subj]]
#		if( sum > 3 ){
			printf p[pos]
			for(subj in s)
				printf ",%s",b[p[pos]][s[subj]]
			printf "\n"
#		}
	}
}
#	This REQUIRES >= gawk 4
