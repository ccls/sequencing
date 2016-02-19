#!/usr/bin/env bash

file=$1
db=$2
log=$file.log
{
echo "Begin $file"
date
hostname
echo

time blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt
#time echo blastn -query $file -db $db -evalue 0.05 -outfmt 0 -out $file.blastn.txt >> $file.log

echo
echo "Done $file"
date
} 1>>$log 2>&1

