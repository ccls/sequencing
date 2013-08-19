#!/bin/sh

echo $0
date

#if [ `uname -n` != 'ec0000' ] ; then
#	echo 
#	echo "You are not on the cluster."
#	echo "You probably don't want to do this."
#	echo
##	exit
#fi
#
#me=`whoami`
#[ ! -z $me ] || me='jwendt'
#echo $me

if [ $# -gt 0 ] ; then
	tmp=`echo $1 | tr -cd '[:digit:]'`
	#	need the x's in case is blank
	if [ "x${tmp}" = "x${1}" ] ; then
		wanted=$1
		shift
#	else
#		wanted=10
	fi
fi
#wanted=${SIMPLE_QUEUE_SCRIPT_SIZE:-10}
[ ! -z $wanted ] || wanted=0

available=`simple_queue.sh size`
[ ! -z $available ] || available=0

echo "Before..."
echo "available ... $available"
echo "wanted    ... $wanted"
echo

if [ $wanted -gt $available ]; then
	echo "More wanted than available."
	wanted=$available
fi

echo "$wanted wanted"
[ ! -z $wanted ] || wanted=0

if [ $wanted -gt 0 ]; then
	now=`date "+%Y%m%d%H%M%S"`
	script_name="queue_script.${now}.sh"
	echo "#!/bin/sh\n\n" > $script_name
	chmod 755 $script_name

	#
	#	seq is different on mac, so wrapping in if
	#
	for i in `seq $wanted` ; do
		echo "simple_queue.sh popping $i"
	#	eval `simple_queue.sh pop`
		echo `simple_queue.sh pop` >> $script_name
	done
fi

echo "After..."
echo "available ... "`simple_queue.sh size`
echo
