#!/bin/sh -x

#	I WANT THE OUTPUT SO KEEP sh -x
##!/usr/bin/env bash


echo
echo
echo "FYI. THIS IS STILL IN DEVELOPMENT!"
echo
echo


QUALITY=${1:-Q20}
INDEX=${2:-hg19}
SOURCE=${3:-/Volumes/box/1000genomes/sequers/1000genomes/untarred/grouping/}

#	This script duplicates the existing references for hg19's Q10 and Q20, but not Q00?

date=`date "+%Y%m%d%H%M%S"`
log_file=`basename $0`.$QUALITY.$INDEX.$date.out

{

find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.$INDEX.$QUALITY.insertion_points.overlappers -depth 2 -exec cat {} \; > $INDEX.$QUALITY.insertion_points.overlappers
find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.both_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points.rc_overlappers -depth 2 -exec cat {} \; > $INDEX.$QUALITY.rc_insertion_points.rc_overlappers

awk '{print $2}' $INDEX.$QUALITY.insertion_points.overlappers | sort -u > $INDEX.$QUALITY.insertion_points.overlappers.sort.uniq
awk '{print $2}' $INDEX.$QUALITY.rc_insertion_points.rc_overlappers | sort -u > $INDEX.$QUALITY.rc_insertion_points.rc_overlappers.sort.uniq


awk '{print $0":F"}' $INDEX.$QUALITY.insertion_points.overlappers.sort.uniq > overlapper_reference.$INDEX.$QUALITY
awk '{print $0":R"}' $INDEX.$QUALITY.rc_insertion_points.rc_overlappers.sort.uniq >> overlapper_reference.$INDEX.$QUALITY

sort overlapper_reference.$INDEX.$QUALITY > overlapper_reference.$INDEX.$QUALITY.sorted

} 1>> $log_file 2>&1

#	not sure what was happening here, but works fine for me on my work laptop?
#	TYPO!  missed an underscore.  empty variable name.  "1>> 2>&1" is bad. 
#	./extract_insertion_points_and_overlappers.sh: line 116: $logfile: ambiguous redirect
#	the {} logfile doesn't seem to work with bash
