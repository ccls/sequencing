#!/usr/bin/env bash

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 ; do
	ls human_genomic-20130527.fa_$n
	bowtie2-build human_genomic-20130527.fa_$n human_genomic_$n
done
