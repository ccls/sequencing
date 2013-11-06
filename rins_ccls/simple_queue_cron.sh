#!/bin/sh

#if [ $# -eq 0 ]; then
#	echo 
#	echo "Will pop and eval commands from simple_queue database"
#	echo "The number will depend on the current content of the slurm queue"
#	echo "as well as the number desired. The default is 200 but this can"
#	echo "be overriden with the SLURM_MAX_QUEUE_SIZE environment variable."
#	echo "This command was designed to be called from a cron job."
#	echo
#	echo "Usage:"
#	echo
#	echo "`basename $0`"
#	echo
#	echo "Example:"
#	echo
#	exit
#fi

echo $0
date

if [ `uname -n` != 'ec0000' ] ; then
	echo 
	echo "You are not on the cluster."
	echo "You probably don't want to do this."
	echo
	exit
fi

me=`whoami`
[ ! -z $me ] || me='jwendt'
echo $me

available=`simple_queue.sh size`
[ ! -z $available ] || available=0
queued=`squeue --noheader --user=$me | wc -l`

#	SLURM_MAX_QUEUE_SIZE MUST BE EXPORTED FROM ENVIRONMENT TO ACTUALLY BE USED
max=${SLURM_MAX_QUEUE_SIZE:-200}
[ ! -z $queued ] || queued=$max

echo "Before..."
echo "available ... $available"
echo "in queue .... $queued"
echo

wanted=`expr $max - $queued`
[ ! -z $wanted ] || wanted=0

if [ $wanted -gt $available ]; then
	echo "More wanted than available."
	wanted=$available
fi

echo "$wanted wanted"
[ ! -z $wanted ] || wanted=0

#
#	Tried with the original while loop,
#	a very similar for loop.
#	Tried removing the eval after removing the trailing " &" from command
#	I just don't get it.  Runs once.  Another fcron instance process is created and stays.
#	Doesn't look like it runs again until the fcron dies?
#
#
#	I still don't get it.  I just doesn't always run.
#	Despite being set to run every 15 minutes, it only runs every 1:15??????  WTF?
#	The cron_test runs every time.
#
#	I have noticed that multiple instances of fcron continue to exist?
#	Kinda like something is causing it to hang.  Commenting these out
#	until I can figure something out. Otherwise, just make it pop more
#	and stop worrying about it.
#
#	jwendt   12466  2452  0 08:55 ?        00:00:00 /usr/sbin/fcron -b
#	jwendt   23180  2452  0 09:14 ?        00:00:00 /usr/sbin/fcron -b
#
#	It is now my understanding that the cron job does not end until everything
#	is completed, even those things run in the background. Until then 
#	/var/log/cron will say "process already running:" and not start it again.
#	In theory, this is nice if you don't want multiple instances running.
#	However, if you don't care, what can you do?
#


for i in `seq $wanted` ; do
	echo "simple_queue.sh popping $i"
	eval `simple_queue.sh pop`
#	cmd=`simple_queue.sh pop | sed 's/ &//'`
#	echo $cmd
#	$cmd &
done

#while [ `simple_queue.sh size` -gt 0 -a `squeue --noheader --user=$me | wc -l` -lt $max ] ; do
#	echo "simple_queue.sh pop"
#	eval `simple_queue.sh pop`
#done

#simple_queue.sh size
#squeue --noheader --user=$me | wc -l

echo "After..."
echo "available ... "`simple_queue.sh size`
echo "in queue .... "`squeue --noheader --user=$me | wc -l`
echo

#	this makes no difference
#exit 0
