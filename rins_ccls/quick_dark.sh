#!/bin/sh -x
#
#	with the -x the commands are sent to STDERR before execution
#
{

#	bowtie output goes to stderr for some reason
#	probably because the SAM file usually goes to stdout
#	so, wrap everything in curly braces and direct both
#	to files.

indexes=/Volumes/cube/working/indexes

bowtie2 -N 1 -q -x $indexes/hg18 \
	-U raw.1.fastq,raw.2.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18.fastq

bowtie2 -N 1 -q -x $indexes/hg19 \
	-U raw_not_hg18.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19.fastq

bowtie2 -N 1 -q -x $indexes/Blast1 \
	-U raw_not_hg18_hg19.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19_Blast1.fastq

bowtie2 -N 1 -q -x $indexes/Blast2 \
	-U raw_not_hg18_hg19_Blast1.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19_Blast1_Blast2.fastq

bowtie2 -N 1 -q -x $indexes/Homo_sapiens.GRCh37.69.cdna.all \
	-U raw_not_hg18_hg19_Blast1_Blast2.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo.fastq

bowtie2 -N 1 -q -x $indexes/nt_human_1 \
	-U raw_not_hg18_hg19_Blast1_Blast2_Homo.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1.fastq 

bowtie2 -N 1 -q -x $indexes/nt_human_2 \
	-U raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1.fastq \
	-S /dev/null --threads 4 \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2.fastq

ifile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2
ofile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2_human_genomic
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 ; do
	ofile=${ofile}_$n
	bowtie2 -N 1 -q -x $indexes/human_genomic_$n \
		-U $ifile.fastq -S /dev/null --threads 4 \
		--un $ofile.fastq
	ifile=$ofile
done

ln -s $ofile.fastq raw_non_human.fastq

#ln -s raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2.fastq raw_non_human.fastq

echo "de novo assembly of single 'unpaired' non-human using Trinity"
Trinity.pl --seqType fq \
	--output trinity_output_single \
	--single raw_non_human.fastq \
	--JM 2G 

cp trinity_output_single/single.fa trinity_input_single.fasta

cp trinity_output_single/Trinity.fasta trinity_non_human_single.fasta

echo "Laning composite fasta file."
bioruby_lane_fasta.rb trinity_input_single.fasta
#		"trinity_input_single_1.fasta".file_check(die_on_failed_file_check)
#		"trinity_input_single_2.fasta".file_check(die_on_failed_file_check)

mv trinity_input_single_1.fasta trinity_input_paired_1.fasta
mv trinity_input_single_2.fasta trinity_input_paired_2.fasta

echo "de novo assembly of re-paired non-human using Trinity"
Trinity.pl --seqType fa \
	--output trinity_output_paired \
	--left  trinity_input_paired_1.fasta \
	--right trinity_input_paired_2.fasta \
	--JM 2G 

cp trinity_output_paired/both.fa trinity_input_paired.fasta

cp trinity_output_paired/Trinity.fasta trinity_non_human_paired.fasta


#
#	Hmm.  Where will the output got when a command redirects its output inside a redirected block?
#
#	Did this in a script and the echo output went to the 'echo' file
#{
#	echo "Echo output" > block_output_test_echo.out 2>&1
#} 1>>block_output_test_block.out 2>&1
#

#blastn -query=trinity_input_single.fasta \
#	-db=$indexes/nt \
#	-evalue 0.05 -outfmt 0 > trinity_input_single_blastn.txt
#
#blastn -query=trinity_non_human_single.fasta \
#	-db=$indexes/nt \
#	-evalue 0.05 -outfmt 0 > trinity_non_human_single_blastn.txt
#
#blastn -query=trinity_input_paired.fasta \
#	-db=$indexes/nt \
#	-evalue 0.05 -outfmt 0 > trinity_input_paired_blastn.txt
#
#blastn -query=trinity_non_human_paired.fasta \
#	-db=$indexes/nt \
#	-evalue 0.05 -outfmt 0 > trinity_non_human_paired_blastn.txt

echo "Finished at ..."
date

#} 1>>quick_dark.out 2>>quick_dark.err
} 1>>quick_dark.out 2>&1
