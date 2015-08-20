#!/bin/bash

#	aws ec2 describe-images --image-ids ami-e77f81a3

image_id="ami-e77f81a3"
instance_type="t2.medium"
volume_size=123


command="aws ec2 run-instances --count 1 \
	--image-id $image_id \
	--instance-type $instance_type \
	--key-name devenv-key \
	--security-groups devenv-sg \
	--instance-initiated-shutdown-behavior terminate \
	--query 'Instances[0].InstanceId' \
	--block-device-mappings '[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":$volume_size}}]'"
echo $command
#instance_id=`$command`
#echo $instance_id
instance_id=`echo $instance_id | sed 's/"//g'`
echo $instance_id

instance_ip=`aws ec2 describe-instances \
	--instance-ids $instance_id \
	--query 'Reservations[0].Instances[0].PublicIpAddress'`
instance_ip=`echo $instance_ip | sed 's/"//g'`
echo $instance_ip


#	start the 1000genomes script in the background.
#	ssh -n -f -i devenv-key.pem ec2-user@$instance_ip "sh -c 'nohup ~/local/bin/aws_1000genomes.sh /dev/null &'"


#	ssh -n -f user@host "sh -c 'cd /whereever; nohup ./whatever > /dev/null 2>&1 &'"
#	ssh askapache 'sh -c "( ( nohup chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'
#	ssh askapache 'nohup sh -c "( ( chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'


