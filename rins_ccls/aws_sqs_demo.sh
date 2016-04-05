#!/usr/bin/env bash


#> aws sqs create-queue --queue-name test123
#{
#    "QueueUrl": "https://us-west-1.queue.amazonaws.com/156714443422/test123"
#}
#> aws sqs get-queue-url --queue-name=test123
#{
#    "QueueUrl": "https://us-west-1.queue.amazonaws.com/156714443422/test123"
#}
#> aws sqs get-queue-url --queue-name=test123  | python -c 'import sys, json; print json.load(sys.stdin)["QueueUrl"]'
#https://us-west-1.queue.amazonaws.com/156714443422/test123


aws sqs send-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123 --message-body "my message body 1"
aws sqs send-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123 --message-body "my message body 2"
aws sqs send-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123 --message-body "my message body 3"
aws sqs send-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123 --message-body "my message body 4"

while true ; do

	message=`aws sqs receive-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123`

	#	once a message is received it is "in flight" and should be deleted when complete.

	echo $message

	if [ -z "$message" ] ; then
		echo "Received message was blank.  I'm done."
		break
	fi
	

	body=`echo $message | python -c 'import sys, json; print json.load(sys.stdin)["Messages"][0]["Body"]'`
	echo $body

	handle=`echo $message | python -c 'import sys, json; print json.load(sys.stdin)["Messages"][0]["ReceiptHandle"]'`
	echo $handle

	#	apparently there is a limit on the number of inflight messages
	#	so delete them as soon as possible.
	aws sqs delete-message --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123 --receipt-handle $handle

	#	if failed, can I unreceive it so as to trigger running at?
	#		or just re-send the message body?

done



#	can only run this once a minute.
#	aws sqs purge-queue --queue-url https://us-west-1.queue.amazonaws.com/156714443422/test123



