#!/bin/sh

if [ $# -eq 0 ]; then
	echo
	echo "maintain FIFO-like sqlite database"
	echo
	echo "Usage:"
	echo "simple_queue.sh pop -> displays and deletes first record in db"
	echo "simple_queue.sh push COMMAND -> adds COMMAND to bottom of db"
	echo "simple_queue.sh size -> display number of records (also count and length)"
	echo "simple_queue.sh list -> display all records"
	echo
	echo "`basename $0`"
	echo
	echo "Example:"
	echo
	exit
fi

#
#	Using -line with sqlite3 doesn't seem to work well within ``.
#	For some reason, it filters out the \n which puts all the output on a 
#	single line which defeats the purpose of the -line parameter.
#
#	Using ~ doesn't work. Use "$HOME" instead.	
#
database_file_name="$HOME/simple_queue.db"
log_file_name="$HOME/simple_queue.log"

#
#	It seems that attempting to read and write to a table when another script is 
#	attempting the same can result in many collisions which result in silent failure.
#	Not entirely silent, but nevertheless.  This is due to the table being locked.
#	Changing the default timeout of 0 to 5000 ms or something will allow the command
#	to wait for the lock to be removed.  This could be done several ways, but must
#	be done again and again with each call.  Apparently it can't be set in the database
#	itself.  It also can't be in the sql command passed to sqlite.  It could either
#	be in a separate init file, or, I just found, as the -cmd value which is sql
#	run before the sql.  
#
#	-init mysqliteinitfile
#
#	-cmd '.timeout 5000' 
#
#	As this should be used in virtually all calls, or at least doesn't hurt
#	those that don't really need it, I'm gonna set a variable here for this.
#
#sqlite="sqlite3 -cmd '.timeout 5000' $database_file_name "
#	using this variable doesn't work?
#sqlite3 -cmd '.timeout 5000' /Users/jakewendt/simple_queue.db
#sqlite3: Error: too many options: "select * from queue"
#Use -help for a list of options.
#sqlite3: Error: too many options: "create table queue(id integer primary key autoincrement, command text)"
#Use -help for a list of options.
#sqlite3: Error: too many options: "select * from queue"
#Use -help for a list of options.
#
#sh-3.2$ sqlite="sqlite3 -cmd '.timeout 5000' $HOME/simple_queue.db "
#sh-3.2$ $sqlite
#Error: unrecognized token: "'.timeout"
#Error: near "/": syntax error
#
#	I tried a number of different quoting, but still won't work?
#

max_delete_retries=5

peek(){
	echo "Peeking ... `date`" >> $log_file_name
	r=`sqlite3 -cmd '.timeout 5000' -line $database_file_name "select * from queue order by id asc limit 1"`
	command=`echo "$r" | grep "^\s*command = " | sed 's/^command = //'`
	echo $command
}

#
#	Sadly, rapid pop calls can result in the multiple returns of the same record.
#
pop(){
	echo "Popping ... `date`" >> $log_file_name
	r=`sqlite3 -cmd '.timeout 5000' -line $database_file_name "select * from queue order by id asc limit 1"`
#	r=`$sqlite -line "select * from queue order by id asc limit 1"`
	#	$r will NOT have the newlines
	# id = 1 command = srun
	#	"$r" WILL have the newlines
	#     id = 1
	#command = srun
#	id=`echo "$r" | grep "^\s*id =" | awk -F= '{print $NF}'`
#	the "^\s*" doesn't work on ec0000?? different grep probably
	id=`echo "$r" | grep " id = " | awk -F= '{print $NF}'`
	if [ "x$id" != "x"  ] ; then

#		command=`echo "$r" | grep "^\s*command = " | awk -F= '{print $NF}'`
#		That don't work if the command has an = in it!!!!!

		command=`echo "$r" | grep "^\s*command = " | sed 's/^command = //'`
		echo $command
		echo $command >> $log_file_name

		#/my/home/jwendt/dna/bin/simple_queue.sh: line 49: 12947 Killed                  sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"
		#	Occassionally, the delete gets killed?  And then is popped and run again.
		#	Not catastrophic, but should be avoided. I would prefer this to not happen, but ...
#		sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"
#		sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"
#		sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"


#
#	would prefer a "do/until" loop so only have the delete command once.
#
#	Actually, could remove the call before the loop.
#	Maybe modify the insides a bit to only echo "failed" if delete_retries > 1 or something
#

		echo "Deleting ..." >> $log_file_name
		sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"
		delete_retries=0
		while [ $delete_retries -lt $max_delete_retries -a \
			`sqlite3 -cmd '.timeout 5000' $database_file_name "select * from queue where id = $id" | wc -l` -gt 0 ]
		do

			#	simple_queue_cron.sh evals the returned value of pop which is
			#	anything that is echoed so don't do this.... unless its to the log file
			echo "Delete failed. Retrying ... $delete_retries" >> $log_file_name

			sqlite3 -cmd '.timeout 5000' $database_file_name "delete from queue where id = $id"
			delete_retries=`expr $delete_retries + 1`
		done
	fi
}

push(){
	#	The $* in a function MUST BE PASSED.  IT IS NOT THE $* from the command line.
	#	Unless, of course, that's what you pass it.
	echo "Pushing ... `date`" >> $log_file_name
	sqlite3 -cmd '.timeout 5000' $database_file_name "insert into queue(command) values('$*')"
	echo "Pushed $*" >> $log_file_name
}

if [ ! -f $database_file_name ] ; then
	sqlite3 -cmd '.timeout 5000' $database_file_name 'create table queue(id integer primary key autoincrement, command text)'
#	$sqlite 'create table queue(id integer primary key autoincrement, command text)'
#else
	#	echo 'found it'
	#
	#	could exist but not be a sqlite db. How to check?
	#
	#	sqlite3 simple_queue.db '.tables'
	#
fi

#	ln is one of the few commands that are atomic.  Using it here to lock the queue
#	to avoid multiple rapid popping to return the same records.
#while ! ln -s ${database_file_name} ${database_file_name}.lock 2>/dev/null ; do :; done
#	add a bit of random sleeping to avoid excessive attempts and continued collisions at linking.
while ! ln -s ${database_file_name} ${database_file_name}.lock 2>/dev/null
	do sleep `echo | gawk '{srand(systime()+PROCINFO["pid"]);print rand()}'` ; done

case "$1" in
	peek )
		shift; peek;;
	pop )
		shift; pop;;
	push )
#		shift; sqlite3 -cmd '.timeout 5000' $database_file_name "insert into queue(command) values('$*')";;
		shift; push $*;;
	size | count | length )
		sqlite3 -cmd '.timeout 5000' $database_file_name "select count(*) from queue" ;;
#	if sqlite does its killing thing, this returns nothing.  Need to fix that.
#	would also like to understand it
	list )
		sqlite3 -cmd '.timeout 5000' $database_file_name "select * from queue";;
esac


host=`hostname`

#	rm is NOT atomic, but mv is
mv ${database_file_name}.lock ${database_file_name}.$$.${host}.deleteme && rm ${database_file_name}.$$.${host}.deleteme
#
#	Actually found an undeleted ${database_file_name}.$$.deleteme
#	Not sure what to do about it.
#	Possible that 2 processes had same PID?  They are running on different computers.
#	Added hostname to filename
#


#while getopts ":db:asdf" clueless
#do
#	echo $clueless
#	echo $OPTARG
#done
#
##	skip to any other command line params
#shift $(($OPTIND - 1 ))
#echo $*
