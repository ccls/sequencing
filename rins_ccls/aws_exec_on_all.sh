#!/bin/bash


function usage(){
	echo
	echo "Usage:"
	echo
	echo "`basename $0` <command>"
	echo
	echo "Example:"
	echo "  `basename $0` 'ls *pid'"
	echo
	exit 1
}
#	Basically, this is TRUE AND DO ...
[ $# -eq 0 ] && usage


for ip in `aws ec2 describe-instances \
	--query 'Reservations[].Instances[].PublicIpAddress' --output text` ; do 

	echo $ip
  ssh -q -n -o UserKnownHostsFile=/dev/null \
		-o StrictHostKeyChecking=no \
		-i devenv-key.pem -l ec2-user $ip $1
#'ls aws_1000genomes.sh.*.pid';

done