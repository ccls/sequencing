This document includes information on installation and use of the PRICE genome assembler.  For information
on running the included sample job, read the accompanying SAMPLE_README.txt document.

For more extensive documentation, download the PRICE documentation HTML pages from the website (or see
the "index.html" file from this download if you downloaded the full set of files).

contact: price@derisilab.ucsf.edu

########## TABLE OF CONTENTS ##########
1. Introduction
2. Installation
3. Basic Use
4. Algorithm Description

########## 1. INTRODUCTION ##########

This is a beta release of the PRICE (Paired Read Iterative Contig Extension) assembler.  It uses paired-end
reads to extend existing contigs.  For metagenomic sequencing projects, sequences can be supplied in a
separate file, possibly a filtered subset of the main datasets, that may derive from an organism/virus of
potential interest and used as seeds to assemble the genome of just that organism/virus using the full
metagenomic dataset.  For single-genome projects, an arbitrary number of reads from the input datasets can
be used to seed assemblies from points scattered across the genome.

PRICE uses paired-end reads to expand contigs.  Sequence data from non-paired-end sources can be used as
initial contigs.

Please be patient with your early beta version of PRICE.  Bugs will certainly be found and corrected in the
coming weeks, documentation will be expanded, and features that are in development will be added.  For the 
latest source code and documentation, check back at derisilab.ucsf.edu for updates.

You can get basic usage information from ./PriceTI with the command-line args -, -h, or --help.

########## 2. INSTALLATION ##########

The .tar file that you just opened should include all of the source code necessary to compile PRICE.  The makefile
will generate a binary PriceTI.  This is the text interface (command-line interface) implementation.  It uses the
gcc compiler and uses all STL components, with the exception of OpenMP, which is used for threading.  OpenMP must
be pre-installed.  For more information on OpenMP, see http://openmp.org/wp/

The following command will compile PRICE:
make

This will create the command-line (i.e. Text-Interface) version of PRICE as a locally-executable binary.  This install
will not find your path and place the executable there, you will need to do that if you want to be able to execute
PRICE from anywhere in your system.  The binary PriceTI has no external dependencies, so it can be copied directly to
other machines with similar architectures.

If you are using a Mac and are having trouble with "make", you probably need to install the Xcode developer tools that
came on a disc with your computer but were not pre-installed.

When installing future versions or updating bug fixes, use the -clean tag before re-compiling:
make clean


########## 3. BASIC USE ##########

After compiling, you can call PRICE from the command line like this:
./PriceTI [args]

Descriptions of the args can be accessed by typing -, -h, or --help as an arg.  Some of the args are explained in greater
detail below.

A sample command:
./PriceTI -fp s_1_1_sequence.txt s_1_2_sequence.txt 300 -picf 100000 s_1_1_sequence.txt 20 3 1 -nc 20 -a 10 -o lane1job.fasta

The output files would be lane1job.cycleN.fasta, where N is every number 1-20.

NOTE: -fp/-fpp can be added multiple times to a command.  If you have multiple lanes' worth of data that you want to
      assemble together (or paired-end data from multiple sequencing platforms), include them all by typing an additional
      -fp flag for each file pair.  The same goes for initial contig files called using -icf/-picf.
NOTE: As in the example above, the same input file can be used repeatedly.  It would not make sense to add the same file
      for multiple -fp/-fpp or multiple -icf/-picf args, but it does make sense to use one of the paired-read files
      that was input using -fp/-fpp as a source of initial contigs using -picf (as shown above).
NOTE: You don't need to write an output file after every cycle.  But I would for now due to the current growth of
      time-per-cycle.
NOTE: Most parts of PRICE are threaded.  -a specifies the number of threads to be used.  Actual thread usage will vary
      during a run.


########## 4. ALGORITHM DESCRIPTION #########

Paired-Read Iterative Contig Extension.  PRICE executes a number of assembly "cycles", each of which has three
steps:
1) mapping of paired-end reads to the edges of existing contigs;
2) assembling the overhanging pairs of mapped reads for each contig edge;
3) performing a meta-assembly of all existing contigs whose aim is to remove redundancy.

### STEP 1 ###
Paired-end reads are provided as file pairs by the -fp and -fpp commands.  The two files must be of equal length,
and the Nth sequence in the first file is taken to be the paired-end of the Nth sequence of the second file.  The
reads are assumed to be from opposite strands of their template molecule and facing one another 5p->3p.  This is the
format of Illumina's s_X_1_sequence.txt and s_X_2_sequence.txt files.  Fasta or fastq (_sequence.txt) files are
accepted.  An additional arg describes the sizes of the amplicons (distance from the 5p edge of one read to the 5p edge
on the opposite strand of its paired end).  It defines an acceptable distance from the contig edge.  Go with the
average amplicon size or a standard deviation or two above that.

  |-----size-------|

5p---->............. 
  .............<----5p

Each pair of files can have a % identity required for mapping to contigs separately specified.  This is useful if
you are mixing datasets from different runs where the quality of the data from those runs is very different.  For
now, this is provided using -fpp instead of -fp.

### STEP 2 ###
Overhanging reads are assembled.  These are reads whose paired-ends mapped to a contig, but they themselves did not
map to any existing contig.  Since they are supposed to extend the contig, that contig itself is also included in
the mini-assembly.  If a paired-end maps to another contig in the dataset, then that contig is included in the 
mini-assembly instead of the read.  This is how contigs from neighboring parts of a genome are brought together.
The amount of sequence identity required to collapse two sequences into one contig is scaled to the number of
sequences in any given mini-assembly.  There are currently arbitrary default values, but they can also be user-adjusted
at the command line.  If you don't like the way your assembly looks, these are the first things to tweak:

Minimum Overlap:
The minimum amount that two sequences need to overlap one another to be combined.
-mol: the global minimum overlap (with fewer sequencs, the min overlap will be decreased, but never below this number).
-tol: the number of sequences at which the global min overlap will be applied (below that number, it will stay constant,
      and above that number, it will be scaled up.  Its increase with respect to the num. sequences will be logarithmic.
      The scaling function is currently not user-controllable.
-mpi: global min % ID
-tpi: threshold for application of min % ID.  Above this number of sequences, the min % ID will be scaled to approach an
      asymptote of 100%.

Initial Contigs:
At the beginning of the assembly, there are no contigs to be extended, so some must be provided.  Generally, I use a small
subset of the input reads.  If you have other, non-paired-end data (like 454 sequence), that is also good, and the longer
the initial sequences are relative to the short read lengths, the faster the assembly will get going.

Initial contigs are taken from files using the -icf and -picf flags.  Using -icf, all of the sequences will be used.
Using -picf, the number N of sequences specified will be used.  They will be taken from throughout the file, not just
the first N sequences.

The initial contigs do not all need to be added at once.  Two numbers are provided to specify the number of steps in which
the initial contigs from that file will be introduced, then the number of cycles that will pass between adding a new batch
of initial contigs.  At this point, I usually add ~1000-5000 at a time for Illumina reads as input (I pick one of the two
read files that I am using, like s_X_1_sequence.txt).

You can also input fasta-formatted contigs from previous assemblies.  But beware: they will be treated as individual reads
unless you use the last arg of -icf/-picf to specify a higher read count.  At the end of each cycle, "contigs" whose max
score count is only 1 are interpreted as leftover individual reads and they are removed.  This is appropriate when they
were just single reads, but is inappropriate when previously-assembled contigs are being used.  In that case, provide a 
number >1 so that the contigs aren't lost.  But don't go too crazy, that will make it hard for them to be modified in
light of conflicting sequence data.  I would go with ~2-5.


### STEP 3 ###
Meta-assembly.  This is really just to remove redundancy, not to put together contigs with small amounts of overlap.
Requirements are still scaled with the number of input sequences just like in the case of the mini-assemblies, but
this is for the sake of computational efficiency.  The commands for scaling are the same as for step 2, but all
uppercase.  Also, note that, since I am going for near-total redundancy, -MOL/-TOL (the overlap args) won't do anything.


