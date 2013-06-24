#!/bin/sh

me=`whoami`

while [ `simple_queue.sh size` -gt 0 -a `squeue --noheader --user=$me | wc -l` -lt 100 ] ; do
	eval `simple_queue.sh pop`
done
