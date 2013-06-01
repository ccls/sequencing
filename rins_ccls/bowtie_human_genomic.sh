#!/bin/sh -x
#
# with the -x the commands are sent to STDERR before execution
#
{
#	quick script to test the removal by human_genomic indexes

ifile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2
ofile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2_human_genomic
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 ; do
	ofile=${ofile}_$n
	bowtie2 -N 1 -q -x /Volumes/cube/working/indexes/human_genomic_$n \
		-U $ifile.fastq -S /dev/null --threads 4 \
		--un $ofile.fastq
	ifile=$ofile
done

} 1>>bowtie_human_genomic.out 2>&1
