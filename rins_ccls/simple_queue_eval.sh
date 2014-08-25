#!/bin/sh

[ -d ~/ec/ ] || mkdir ~/ec/
[ -d ~/ec/logs/ ] || mkdir ~/ec/logs/
[ -d ~/ec/pids/ ] || mkdir ~/ec/pids/


#       leading with the ": " stops execution
#       just ${BLASTDB:"/Volumes/cube/working/indexes"}
#       would try to execute the result.  I just want the OR/EQUALS feature
: ${SLURMD_NODENAME:="n0000"}
: ${SLURM_JOBID:="000000"}
: ${SLURM_TASK_PID:="0000"}

touch ~/ec/pids/$SLURMD_NODENAME.$SLURM_JOBID.$SLURM_TASK_PID

#	SLURM_NODELIST=n0006
#	SLURM_JOBID=231351
#	SLURM_JOB_NAME=slurm_env.sh
#	SLURMD_NODENAME=n0006
#	SLURM_TASK_PID=13731

{
	echo "Starting at ..."
	date

	while [ -f ~/ec/pids/$SLURMD_NODENAME.$SLURM_JOBID.$SLURM_TASK_PID ]; do

		#
		#	FYI. Counting can kill as it isn't dealt with if sqlite has a hissy fit.
		#		Seems that if sqlite does its killing thing, returns blank.
		#
		available=`simple_queue.sh size`
		available_status=$?
		[ ! -z $available ] || available=0

		if [ $available -gt 0 ]; then

			cmd=`simple_queue.sh pop`
			echo "- - - - - - - - - -"
			echo "Popped ..."
			echo $cmd
			#
			#	remove everything before blastn
			#		srun --share --job-name=003042_viral_genomic 
			#
			#	AND remove the trailing &!!!!!
			#
			#	Queue should no longer contain these.
#			cmd=`echo $cmd | sed 's/^srun.*blastn -q/blastn -q/' | sed 's/ &$//'`
			echo "Executing ..."
			echo $cmd
			date
			output=`$cmd`
			status=$?
			if [ $status -ne 0 ]; then
				echo $cmd >> ~/ec/logs/$SLURMD_NODENAME.$SLURM_JOBID.$SLURM_TASK_PID.FAILED
			fi
			echo "Completed. :${output}:"	#	"output" will always be blank
			date
			echo "- - - - - - - - - -"
		else
			echo "Available :${available}: not greater than zero.  Status:${available_status}:. Commiting suicide.  Goodbye cruel world!"
			rm ~/ec/pids/$SLURMD_NODENAME.$SLURM_JOBID.$SLURM_TASK_PID
		fi

	done

	echo "PID file has been removed"
	echo "Ending at ..."
	date

} 1>>~/ec/logs/$SLURMD_NODENAME.$SLURM_JOBID.$SLURM_TASK_PID.log 2>&1

exit;

#	Run this 4 times.  14 seems to be the max callable and 4*14=56.  Conveniently.
#	srun --share --partition=ccls --ntasks=14 simple_queue_eval.sh &

#	There are actually 56 run units available on ccls (not 52)
#	n0003 - 16
#	n0004 - 16
#	n0005 - 16
#	n0006 - 4
#	n0007 - 4
