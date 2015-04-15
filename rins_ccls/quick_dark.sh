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


if [ $# -ne 2 ]; then
	echo
	echo "Usage:"
	echo
	echo "`basename $0` raw.1.fastq raw.2.fastq"
	echo
	echo "Example:"
	echo
	exit
fi

{
echo "Starting at ..."
date

#if [ $# -eq 2 ] ; then
	ln -s $1 raw.1.fastq
	ln -s $2 raw.2.fastq
#fi
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

bowtie2="bowtie2 -N 1 -q -S /dev/null --threads 4"

#
#	Try to super simplify the bowtie stuff like so ....
#
#	nt_human_1 nt_human_2,2 nt_human_3,3
dbs="hg18 hg19 Blast1 Blast2,2 Homo_sapiens.GRCh37.69.cdna.all,Homo
	human_genomic_01
	human_genomic_02,02
	human_genomic_03,03
	human_genomic_04,04
	human_genomic_05,05
	human_genomic_06,06
	human_genomic_07,07
	human_genomic_08,08
	human_genomic_09,09
	human_genomic_10,10
	human_genomic_11,11
	human_genomic_12,12
	human_genomic_13,13"

ifile='raw.1.fastq,raw.2'	#	yes, this works
ofile='raw_not'
for db in $dbs; do
	#	If contains comma (nt_human_2,2), split it.  If not (hg18), db and suffix will be the same.
	suffix=${db##*,}	#	2
	db=${db%%,*}			#	nt_human_2
	ofile=${ofile}_$suffix
	$bowtie2 -x $db -U $ifile.fastq --un $ofile.fastq

	status=$?
	if [ $status -ne 0 ] ; then
		date
		echo "bowtie failed with $status"
		exit $status
	fi

	ifile=$ofile
	date
done


#ln -s $ofile.fastq raw_non_human.fastq
mv $ofile.fastq raw_non_human.fastq


#	to speed things onto the cluster, convert input fastq to fasta before trinity.
#
#	Looks like the trinity fastq converter changes the sequence name whilst the fastx command does not.
#	Not sure if this will be an issue, but noting it. (had to modify my script bioruby_lane_fasta.rb)
#
#	==> raw_non_human.fasta <==
#	>HWI-700460R:370:C38TNACXX:8:1101:12930:2231 1:N:0:CGATGT
#	==> trinity_input_single.fasta <==
#	>HWI-700460R:370:C38TNACXX:8:1101:12930:2231/1
#
#	The -n means keep sequences with unknown (N) nucleotides. (Do we want that?)
#	-Q33 is UNDOCUMENTED AND NEEDED for our fastq files.

echo
echo "Converting FASTQ raw_non_human.fastq to FASTA raw_non_human.fasta"
fastq_to_fasta -Q33 -n -i raw_non_human.fastq -o trinity_input_single.presed.fasta
status=$?

#	sometimes bowtie mucks up the file during processing
#	I'm guessing that this is the last run as bowtie did not, but has in the past, complain
#fastq_to_fasta: Error: invalid quality score data on line 58491284 (quality_tok = "@C@FFFDFHHGHGJIJJIIJHGJJIJJIJIJJJIGIIIIIIJIJEHHFHFDFFFFEEEEECDDBCD@DACCDBBDDDD@>BDDCCCD@HWI-700460R:370:C38TNACXX:8:1308:8501:2767 1:N:0:ACAGTG"

#	for some reason, the last 13 chars "CDABDDCD@<>9<" of the previous quality score were tossed?
#	This really makes me question stuff.  What else is wrong, but not syntactically correct?
#	Is it bowtie?  The filesystem?  The disk?  The network?
#
#	Do not delete these files until the very end.
#	Checking if bowtie will complain on the corruption.
#	fortunately bowtie will fail on this so certain the last bowtie caused this.
#	Error: Encountered one or more spaces while parsing the quality string for read HWI-700460R:370:C38TNACXX:8:1308:7929:2821 1:N:0:ACAGTG.  If this is a FASTQ file with integer (non-ASCII-encoded) qualities, try re-running with the --integer-quals option.
#	libc++abi.dylib: terminating with uncaught exception of type int
#	bowtie2-align died with signal 6 (ABRT) 
#

#@HWI-700460R:370:C38TNACXX:8:1308:7929:2821 1:N:0:ACAGTG
#GGAGGGAGGAAGACGAACGGAAGGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA
#+
#@C@FFFDFHHGHGJIJJIIJHGJJIJJIJIJJJIGIIIIIIJIJEHHFHFDFFFFEEEEECDDBCD@DACCDBBDDDD@>BDDCCCD@HWI-700460R:370:C38TNACXX:8:1308:8501:2767 1:N:0:ACAGT        G
#ACGAACGGAAGGACGGACGGCGCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAA
#+
#@@<DDDDADCFHB;E6F::?8@<FBGE)=85?BBBD>@CBB(55=>>C:@C3@@A@@CCBBA:@(+>@<8?>3::99<@B>4:>A::1:A@#########
#
#[jwendt@n0 MeSF001_SF3]$ grep -A 3  "@HWI-700460R:370:C38TNACXX:8:1308:7929:2821" raw.1.fastq 
#
#@HWI-700460R:370:C38TNACXX:8:1308:7929:2821 1:N:0:ACAGTG
#GGAGGGAGGAAGACGAACGGAAGGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA
#+
#@C@FFFDFHHGHGJIJJIIJHGJJIJJIJIJJJIGIIIIIIJIJEHHFHFDFFFFEEEEECDDBCD@DACCDBBDDDD@>BDDCCCDCDABDDCD@<>9<


#	so, oddly enough, fastq_to_fasta failed again on line 205370349
#	The missing characters from above ended up about 150000000 lines later?????
#	But not in the original file, just the one after sed fixed the first problem????
#	hidden control characters?
#	MY FAULT!  I didn't include the lane in the sed so it did it twice as it was told.  My bad.
#	Just the first problem.

if [ $status -ne 0 ] ; then
	date
	echo "fastq_to_fasta failed with $status"
	exit $status
fi

#	The output will be re-paired and run through Trinity which seems to REQUIRE
#	that the sequences names be laned with /1 or /2 and NOT like ... 1:N:0:CGGAAT
#warning, ignoring read: HS2:360:D1NTUACXX:8:1101:12056:83547 since cannot decipher if /1 or /2 of a pair.
#	So we can use Trinity's converter, another converter or convert the output here?
#	The RINS script doesn't do it.
#	Trinity's fastool is buried inside it and not in the path.
#tail trinity_input_single.fasta | sed 's;\(>.*\) \([12]\):.*$;\1/\2;'
#jakewendt@fxdgroup-169-229-196-225 : TTV7 929> tail trinity_input_single.fasta | sed 's;\(>.*\) \([12]\):.*$;\1/\2;'
#>HS2:360:D1NTUACXX:8:2308:19573:200502 2:Y:0:CGGAAT
#>HS2:360:D1NTUACXX:8:2308:19573:200502/2

#	Also, fastx_collapser TOTALLY mucks up the read counts if the sequence name has a - in it.
#	>HWI-700460R:370:C38TNACXX:8:1101:8490:2221 1:N:0:ATNACG
#	there should be no - anywhere else in the fasta file so just replace all.
#
#	FYI, I don't think this is actually a muck up.  I think that it is trying to preserve
#	a previous counts by adding the counts in the sequence names by the number following the first
#	dash in the sequence name, if there is one.  If too large, the number loops over and becomes 
#	negative.  When really large, It loops around a couple times.
#	It would be nice if there was some documentation about this.
#















#	TAG EACH READ WITH THE SAMPLE NAME?
#	Names can include dashes as well as other chars so replace them with underscores as well.















sed 's;\(>.*\) \([12]\):.*$;\1/\2;' trinity_input_single.presed.fasta | sed 's/-/_/g' > trinity_input_single.fasta
status=$?
if [ $status -ne 0 ] ; then
	date
	echo "sed failed with $status"
	exit $status
fi








echo
echo "Removing duplicate reads from fasta files to speed up blasting."

#bioruby_extract_uniq_sequences_from_fasta.rb trinity_input_single.fasta
#	=> trinity_input_single.uniq.fasta
#	using fastx_collapser instead.  Better on memory and faster.
#	Also adds read count
#	http://hannonlab.cshl.edu/fastx_toolkit/
#	This will completely rename the reads so will lose lane info etc.
fastx_collapser -i trinity_input_single.fasta -o trinity_input_single.uniq.fasta
status=$?
if [ $status -ne 0 ] ; then
	date
	echo "fastx_collapser failed with $status"
	exit $status
fi




echo
echo "Splitting input fasta file into 20000 read fasta files and" \
	"queueing for blasting to viral genomic"
date
ec_fasta_split_and_blast.sh --max_reads 20000 \
	--dbs viral_genomic \
	--options "-task blastn" \
	trinity_input_single.uniq.fasta

#	This sleep is used to ensure that the directory created above
#	does not have the same timestamp as the one below.
sleep 2

echo
echo "Splitting input fasta file into 2000 read fasta files and" \
	"queueing for blasting to nt"
date
ec_fasta_split_and_blast.sh --max_reads 2000 \
	trinity_input_single.uniq.fasta








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
#	20141009 - added --run_as_paired.  Not sure how this will affect anything.
#             (we primarily use the paired output anyway)
#	20150303 - change --JM to --max_memory
#

date
echo
echo "de novo assembly of single 'unpaired' non-human using Trinity"
Trinity --seqType fa --bflyHeapSpaceMax 5G --max_memory 2G \
	--run_as_paired \
	--min_contig_length 100 \
	--single trinity_input_single.fasta \
	--output trinity_output_single
date

#cp trinity_output_single/single.fa trinity_input_single.fasta

cp trinity_output_single/Trinity.fasta trinity_non_human_single.fasta


echo
echo "Splitting output fasta file into 20000 read fasta files and" \
	"queueing for blasting to viral genomic"
date
ec_fasta_split_and_blast.sh --max_reads 20000 \
	--dbs viral_genomic \
	--options "-task blastn" \
	trinity_non_human_single.fasta

#	This sleep is used to ensure that the directory created above
#	does not have the same timestamp as the one below.
sleep 2

echo
echo "Splitting output fasta file into 2000 read fasta files and" \
	"queueing for blasting to nt"
date
ec_fasta_split_and_blast.sh --max_reads 2000 \
	trinity_non_human_single.fasta



echo
echo "Laning composite fasta file."
bioruby_lane_fasta.rb trinity_input_single.fasta
#	=> trinity_input_single_R1.fasta, trinity_input_single_R2.fasta

mv trinity_input_single_R1.fasta trinity_input_paired_R1.fasta
mv trinity_input_single_R2.fasta trinity_input_paired_R2.fasta

echo
echo "de novo assembly of re-paired non-human using Trinity"
Trinity --seqType fa --bflyHeapSpaceMax 5G --max_memory 2G \
	--min_contig_length 100 \
	--left  trinity_input_paired_R1.fasta \
	--right trinity_input_paired_R2.fasta \
	--output trinity_output_paired
date

#
#	We are no longer keeping trinity_input_paired related files (subset of trinity_input_single)
#
#	cp trinity_output_paired/both.fa trinity_input_paired.fasta

cp trinity_output_paired/Trinity.fasta trinity_non_human_paired.fasta


#
#	We are no longer keeping trinity_input_paired related files (subset of trinity_input_single)
#
#bioruby_extract_uniq_sequences_from_fasta.rb trinity_input_paired.fasta
#	=> trinity_input_paired.uniq.fasta


#	even 10000 processes quite fast
echo
echo "Splitting output fasta file into 20000 read fasta files and" \
	"queueing for blasting to viral genomic"
date
ec_fasta_split_and_blast.sh --max_reads 20000 \
	--dbs viral_genomic \
	--options "-task blastn" \
	trinity_non_human_paired.fasta

#	This sleep is used to ensure that the directory created above
#	does not have the same timestamp as the one below.
sleep 2

#	Defaults are -m 1000 and --dbs nt
echo
echo "Splitting output fasta file into 2000 read fasta files and" \
	"queueing for blasting to nt"
date
ec_fasta_split_and_blast.sh --max_reads 2000 \
	trinity_non_human_paired.fasta


#	Defaults are -m 1000 and --dbs nt
echo
echo "Splitting output fasta file into 10000 read fasta files" \
	"and queueing for tblastx'ing to viral_genomic"
date
ec_fasta_split_and_blast.sh --command tblastx --max_reads 10000 \
	trinity_non_human_paired.fasta


#
#	We are no longer keeping trinity_input_paired related files (subset of trinity_input_single)
#
#ec_fasta_split_and_blast.sh trinity_input_paired.uniq.fasta


echo
echo "Finished at ..."
date

#} 1>>quick_dark.out 2>>quick_dark.err
# With the -x the commands are sent to STDERR before execution
#	so send both STDOUT and STDERR to same file.
} 1>>quick_dark.out 2>&1
