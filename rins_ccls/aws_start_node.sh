#!/bin/bash

#	aws ec2 describe-images --image-ids ami-4335cc07

image_id="ami-4335cc07"
instance_type="t2.medium"
#volume_size=123


instance_id=`aws ec2 run-instances --count 1 \
	--image-id $image_id \
	--instance-type $instance_type \
	--key-name devenv-key \
	--security-groups devenv-sg \
	--instance-initiated-shutdown-behavior terminate \
	--query 'Instances[0].InstanceId'`
#	--block-device-mappings '[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":$volume_size}}]'"
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


#	start the 1000genomes script in the background.
#	if path is set in bashrc (instead of bash_profile), don't need path
#	ssh -n -f -i devenv-key.pem ec2-user@$instance_ip 'sh -c "( ( nohup aws_1000genomes.sh &>/dev/null ) & )"'


#	ssh -n -f user@host "sh -c 'cd /whereever; nohup ./whatever > /dev/null 2>&1 &'"
#	ssh askapache 'sh -c "( ( nohup chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'
#	ssh askapache 'nohup sh -c "( ( chown -R ask:ask /www/askapache.com &>/dev/null ) & )"'


