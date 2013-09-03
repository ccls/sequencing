#!/bin/sh -x
#
#	with the -x the commands are sent to STDERR before execution
#
#	bowtie output goes to stderr for some reason
#	probably because the SAM file usually goes to stdout
#	so, wrap everything in curly braces and direct both
#	to files.
#
#	Explicit redirection within the block will override this
#


if [ $# -eq 0 ]; then
	echo
	echo "Usage:"
	echo
	echo "`basename $0` raw.1.fastq and raw.2.fastq"
	echo
	echo "Example:"
	echo
	exit
fi

{
echo "Starting at ..."
date

#if [ $# -eq 1 ] ; then
#	echo "1 arg."
#elif [ $# -eq 2 ] ; then
#	echo "2 args."
#elif [ $# -gt 2 ] ; then
#	echo "More than 2 args."
#else
#	echo "Zero args."
#fi

if [ $# -eq 2 ] ; then
	ln -s $1 raw.1.fastq
	ln -s $2 raw.2.fastq
fi
#	else ...
#	I expect the existance of raw.1.fastq and raw.2.fastq
#

#indexes=/Volumes/cube/working/indexes
#	leading with the ": " stops execution
#	just ${BOWTIE2_INDEXES:"/Volumes/cube/working/indexes"}
#	would try to execute the result.  I just want the OR/EQUALS feature
: ${BOWTIE2_INDEXES:="/Volumes/cube/working/indexes"}
: ${BLASTDB:="/Volumes/cube/working/indexes"}

#	they MUST be exported, apparently, to be picked up by bowtie2
export BOWTIE2_INDEXES
export BLASTDB

#bowtie2
#-x <bt2-idx> The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.

#	blastn
#$BLASTDB - This is the variable which points to the Blast Database. This directory should contain the databases that you would want to search. BLAST by default checks this location and the current working directory for the presence of the databases. This variable is set during login by system login scripts , and may be changed by the user to point to her preferred location in her startup scripts. 

bowtie2="bowtie2 -N 1 -q -S /dev/null --threads 4 "

$bowtie2 -x hg18 \
	-U raw.1.fastq,raw.2.fastq \
	--un raw_not_hg18.fastq

$bowtie2 -x hg19 \
	-U raw_not_hg18.fastq \
	--un raw_not_hg18_hg19.fastq

$bowtie2 -x Blast1 \
	-U raw_not_hg18_hg19.fastq \
	--un raw_not_hg18_hg19_Blast1.fastq

$bowtie2 -x Blast2 \
	-U raw_not_hg18_hg19_Blast1.fastq \
	--un raw_not_hg18_hg19_Blast1_Blast2.fastq

$bowtie2 -x Homo_sapiens.GRCh37.69.cdna.all \
	-U raw_not_hg18_hg19_Blast1_Blast2.fastq \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo.fastq



$bowtie2 -x nt_human_1 \
	-U raw_not_hg18_hg19_Blast1_Blast2_Homo.fastq \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1.fastq 

$bowtie2 -x nt_human_2 \
	-U raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1.fastq \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2.fastq

$bowtie2 -x nt_human_3 \
	-U raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2.fastq \
	--un raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2_3.fastq



ifile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2_3
ofile=raw_not_hg18_hg19_Blast1_Blast2_Homo_nt_human_1_2_3_human_genomic
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 ; do
	ofile=${ofile}_$n
	$bowtie2 -x human_genomic_$n -U $ifile.fastq --un $ofile.fastq
	ifile=$ofile
done


ln -s $ofile.fastq raw_non_human.fastq


#
#	On genepi/n0, had some failure
#
#CMD: /my/home/jwendt/dna/trinityrnaseq_r2013-02-25/trinity-plugins/parafly/bin/ParaFly -c /my/home/jwendt/dna/output/fallon_715723_filtered/trinity_output_single/chrysalis/butterfly_commands.adj -shuffle -CPU 2 -failed_cmds failed_butterfly_commands.20633.txt -v 
#Number of  Commands: 3057
#Error occurred during initialization of VM
#Could not reserve enough space for object heap
#Error: Could not create the Java Virtual Machine.
#Error: A fatal exception has occurred. Program will exit.
#
#	Adding --bflyHeapSpaceMax 5G to Trinity.pl call seems to fix
#
#	adding --min_contig_length 100 in attempt to get ALL input reads in the output
#

echo "de novo assembly of single 'unpaired' non-human using Trinity"
Trinity.pl --seqType fq --bflyHeapSpaceMax 5G --JM 2G \
	--min_contig_length 100 \
	--single raw_non_human.fastq \
	--output trinity_output_single

cp trinity_output_single/single.fa trinity_input_single.fasta

cp trinity_output_single/Trinity.fasta trinity_non_human_single.fasta

echo "Laning composite fasta file."
bioruby_lane_fasta.rb trinity_input_single.fasta

mv trinity_input_single_1.fasta trinity_input_paired_1.fasta
mv trinity_input_single_2.fasta trinity_input_paired_2.fasta

echo "de novo assembly of re-paired non-human using Trinity"
Trinity.pl --seqType fa --bflyHeapSpaceMax 5G --JM 2G \
	--min_contig_length 100 \
	--left  trinity_input_paired_1.fasta \
	--right trinity_input_paired_2.fasta \
	--output trinity_output_paired

cp trinity_output_paired/both.fa trinity_input_paired.fasta

cp trinity_output_paired/Trinity.fasta trinity_non_human_paired.fasta


#
#	This is where I would like to scp a couple fasta files to the cluster
#		and then blast 'em with ...
#	ssh -f genepi ssh -f ec0000 blastn ....
#	... BUT ...
#



#
#	Hmm.  Where will the output got when a command redirects its output inside a redirected block?
#
#	Did this in a script and the echo output went to the 'echo' file
#{
#	echo "Echo output" > block_output_test_echo.out 2>&1
#} 1>>block_output_test_block.out 2>&1
#

#blastn -query=trinity_input_single.fasta \
#	-db=nt \
#	-evalue 0.05 -outfmt 0 > trinity_input_single_blastn.txt
#
#blastn -query=trinity_non_human_single.fasta \
#	-db=nt \
#	-evalue 0.05 -outfmt 0 > trinity_non_human_single_blastn.txt
#
#blastn -query=trinity_input_paired.fasta \
#	-db=nt \
#	-evalue 0.05 -outfmt 0 > trinity_input_paired_blastn.txt
#
#blastn -query=trinity_non_human_paired.fasta \
#	-db=nt \
#	-evalue 0.05 -outfmt 0 > trinity_non_human_paired_blastn.txt

echo "Finished at ..."
date

#} 1>>quick_dark.out 2>>quick_dark.err
} 1>>quick_dark.out 2>&1
