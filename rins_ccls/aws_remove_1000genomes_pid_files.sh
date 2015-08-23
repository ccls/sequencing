#!/bin/bash

for ip in `aws ec2 describe-instances \
	--query 'Reservations[].Instances[].PublicIpAddress' --output text` ; do 
  ssh -q -n -o UserKnownHostsFile=/dev/null \
		-o StrictHostKeyChecking=no \
		-i devenv-key.pem -l ec2-user $ip 'uname -a';
done
