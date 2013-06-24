#!/bin/sh
#
#	Using -line with sqlite3 doesn't seem to work well within ``.
#	For some reason, it filters out the \n which puts all the output on a 
#	single line which defeats the purpose of the -line parameter.
#
#	Using ~ doesn't work. Use "$HOME" instead.	
#
database_file_name="$HOME/simple_queue.db"

pop(){
	r=`sqlite3 -line $database_file_name "select * from queue order by id asc limit 1"`
	#	$r will NOT have the newlines
	# id = 1 command = srun
	#	"$r" WILL have the newlines
	#     id = 1
	#command = srun
	id=`echo "$r" | grep "^\s*id =" | awk -F= '{print $NF}'`
	if [ "x$id" != "x"  ] ; then
		command=`echo "$r" | grep "^\s*command =" | awk -F= '{print $NF}'`
		echo $command
		sqlite3 $database_file_name "delete from queue where id = $id"
	fi
}


if [ ! -f $database_file_name ] ; then
	sqlite3 $database_file_name 'create table queue(id integer primary key autoincrement, command text)'
#else
	#	echo 'found it'
	#
	#	could exist but not be a sqlite db. How to check?
	#
	#	sqlite3 simple_queue.db '.tables'
	#
fi

case "$1" in
	pop )
		shift; pop;;
	push )
		shift; sqlite3 $database_file_name "insert into queue(command) values('$*')";;
	size | count | length )
		sqlite3 $database_file_name "select count(*) from queue" ;;
	* )
		sqlite3 $database_file_name "select * from queue"
esac


#while getopts ":db:asdf" clueless
#do
#	echo $clueless
#	echo $OPTARG
#done
#
##	skip to any other command line params
#shift $(($OPTIND - 1 ))
#echo $*
