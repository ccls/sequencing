#!/bin/sh -x

QUALITY=${1:-Q20}
INDEX=${2:-hg19}
SOURCE=${3:-/Volumes/box/1000genomes/sequers/1000genomes/untarred/grouping/}

date=`date "+%Y%m%d%H%M%S"`
log_file=`basename $0`.$QUALITY.$INDEX.$date.out

{

find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.$INDEX.$QUALITY.insertion_points -depth 2 -exec cat {} \; > pre_ltr.bowtie2.$INDEX.$QUALITY.insertion_points
find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.pre_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points -depth 2 -exec cat {} \; > pre_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points
find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.post_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points -depth 2 -exec cat {} \; > post_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points
find $SOURCE -name \*.bowtie2.herv_k113_ltr_ends.__very_sensitive_local.aligned.post_ltr.bowtie2.$INDEX.$QUALITY.insertion_points -depth 2 -exec cat {} \; > post_ltr.bowtie2.$INDEX.$QUALITY.insertion_points

sort -u pre_ltr.bowtie2.$INDEX.$QUALITY.insertion_points > pre_ltr.bowtie2.$INDEX.$QUALITY.insertion_points.sort.uniq
sort -u pre_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points > pre_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points.sort.uniq
sort -u post_ltr.bowtie2.$INDEX.$QUALITY.insertion_points > post_ltr.bowtie2.$INDEX.$QUALITY.insertion_points.sort.uniq
sort -u post_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points > post_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points.sort.uniq

positions_within_10bp.sh {post,pre}_ltr.bowtie2.$INDEX.$QUALITY.insertion_points.sort.uniq | sort -u > overlapper_reference.$INDEX.$QUALITY.F
positions_within_10bp.sh {post,pre}_ltr.bowtie2.$INDEX.$QUALITY.rc_insertion_points.sort.uniq | sort -u > overlapper_reference.$INDEX.$QUALITY.R

awk '{print $0":F"}' overlapper_reference.$INDEX.$QUALITY.F > overlapper_reference.$INDEX.$QUALITY
awk '{print $0":R"}' overlapper_reference.$INDEX.$QUALITY.R >> overlapper_reference.$INDEX.$QUALITY


} 1>>$log_file 2>&1
#	the {} logfile doesn't seem to work with bash
