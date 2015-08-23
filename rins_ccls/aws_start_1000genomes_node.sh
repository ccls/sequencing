#!/bin/bash

#	aws ec2 describe-images --image-ids ami-4335cc07

image_id="ami-b14eb7f5"
#instance_type="t2.micro"
instance_type="t2.medium"
#volume_size=10


instance_id=`aws ec2 run-instances --count 1 \
	--image-id $image_id \
	--instance-type $instance_type \
	--key-name devenv-key \
	--security-groups devenv-sg \
	--instance-initiated-shutdown-behavior terminate \
	--user-data file://aws_start_1000genomes_processing.sh \
	--query 'Instances[0].InstanceId'`
#	--block-device-mappings '[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":$volume_size}}]'

echo $instance_id
instance_id=`echo $instance_id | sed 's/"//g'`
echo $instance_id

echo "Waiting a moment to get an ip address ..."
sleep 5

instance_ip=`aws ec2 describe-instances \
	--instance-ids $instance_id \
	--query 'Reservations[0].Instances[0].PublicIpAddress'`
instance_ip=`echo $instance_ip | sed 's/"//g'`
echo $instance_ip

echo ssh -i devenv-key.pem ec2-user@$instance_ip

#echo "Waiting 120 seconds to for instance to start ..."
#sleep 120
#
##	start the 1000genomes script in the background.
##	if path is set in bashrc (instead of bash_profile), don't need path
##	ssh -n -f -i devenv-key.pem ec2-user@$instance_ip 'sh -c "( ( nohup aws_1000genomes.sh &>/dev/null ) & )"'
##		won't sudo shutdown script because of some tty issue.
##		
##		sudo the whole script?
##	ssh -n -f -i devenv-key.pem ec2-user@$instance_ip 'sh -c "( ( sudo nohup aws_1000genomes.sh &>/dev/null ) & )"'
#
#
#echo "Here ... we ... go!"
#echo ssh -n -f -i devenv-key.pem ec2-user@$instance_ip 'sh -c "( ( nohup aws_1000genomes.sh --shutdown &> aws_1000genomes.sh.log ) & )"'
#ssh -n -f -i devenv-key.pem ec2-user@$instance_ip 'sh -c "( ( nohup aws_1000genomes.sh --shutdown &> aws_1000genomes.sh.log ) & )"'
#
#
##	ssh -n -f user@host "sh -c 'cd /whereever; nohup ./whatever > /dev/null 2>&1 &'"
##	ssh askapache 'sh -c "( ( nohup chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'
##	ssh askapache 'nohup sh -c "( ( chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'
#
#
##	sudo visudo 
##		to comment out the following lines ...
##	Defaults    requiretty
##	Defaults   !visiblepw
#
