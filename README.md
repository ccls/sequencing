# CCLS Sequencing

This is simply a collection of notes and scripts.
And now, other developer's source code.


[RINS Package](http://khavarilab.stanford.edu/resources.html)

Eventually, I'd prefer to understand exactly what RINS and the
collection of apps that it uses does and then replicate in ruby
with testing.  RINS includes some examples, but none have
worked for me, so far.

	Parameters have changed.  
	Databases are corrupt.
	Executables are seg faulting.
	Processing "chopped" files finds nothing.
	Parafly java crashes...
	 	( removed --compatible_path_extension)
	blastn output changed?
		modified blastn_cleanup.pl and write_results.pl




I have begun including source code and making modifications to allow
for the complete install of Blat, Blast, Bowtie, Bowtie2, Trinity and RINS
from this repository.  It should be clearly noted that the aforementioned
software packages ARE NOT MINE and belong to those that wrote them.

I created the root Makefile to control the making and installing of these 
packages on MY MACHINE (Mac Pro 10.6.8) and so far this works.
Used GCC 421



[RINS Package](http://khavarilab.stanford.edu/resources.html) 
Not including the virus index.

[RINS Data (not included)](https://s3.amazonaws.com/changseq/kqu/rins/rins.tar.gz)

[BLAST 2.2.27](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
( 'make'ing this generates about 4000 additional files.  Odd. )
 
[Blat 34](http://users.soe.ucsc.edu/~kent/src/) 
THIS IS NOT THE LATEST VERSION.

[Blat 35](http://users.soe.ucsc.edu/~kent/src/) 

[Bowtie 0.12.8 and Bowtie 2.0.0 beta 7](http://bowtie-bio.sourceforge.net/)

[Trinity 2012-06-08](http://trinityrnaseq.sourceforge.net) 
Not including the ~50MB of sample data

[Trinity 2012-10-05](http://trinityrnaseq.sourceforge.net) 
Not including the ~50MB of sample data

[BWA 0.6.2](http://sourceforge.net/projects/bio-bwa/files/)

[PRICE 20120527](http://derisilab.ucsf.edu/software/price/)

[Velvet](http://www.ebi.ac.uk/~zerbino/velvet/)

[VICUNA](http://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/vicuna)

[MIRA](http://www.chevreux.org/projects_mira.html)


[BioRuby](http://bioruby.org/) - not yet used for anything



## Indexes Used

Manually downloaded, but then for some reason couldn't actually use?
ftp://ftp.ncbi.nlm.nih.gov/blast/db/
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz


So used update\_blastdb.pl to get them

	update_blastdb.pl --decompress nt
	update_blastdb.pl --decompress nr

This will download them to whereever you are so make sure you have ~50GB of space!


	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
	twoBitToFa hg19.2bit hg19.fa
	#       generate ebwt files
	bowtie-build hg19.fa hg19

The big rins.tar.gz includes hg18 and virus


## General modifications

Changed all perl scripts she-bang to use environment selected perl.
Blat doesn't have any perl, so this was just Bowtie, Bowtie2, Blast and Trinity.

	From ...

		#!/usr/bin/perl

	... to ...

		#!/usr/bin/env perl


## Included BLAT 34 modifications required for me to compile

	blatSrc34/jkOwnLib/gfPcrLib.c
	  errAbort(buf);
	    ... changed to ...
	  errAbort("%s", buf);

	blatSrc34/lib/errCatch.c
	  warn(errCatch->message->string);
	    ... changed to ...
	  warn("%s",errCatch->message->string);

	blatSrc34/lib/htmlPage.c
	  warn(errCatch->message->string);
	    ... changed to ...
	  warn("%s", errCatch->message->string);



## Trinity 2012-06-08 modifications ( scripts to be copied in via Makefile )

SAM\_filter\_out\_unmapped\_reads.pl was modified due to some test
data triggering a divide by zero error.

Add an empty target to a Makefile.in.
I doesn't seem needed on my Mac Pro, just my MacBook Pro? 
Different version of make? XCode?

	> trinityrnaseq_r2012-06-08/trinity-plugins/jellyfish-1.1.5/Makefile.in 

		doc/jellyfish.man:
		#      this is just an empty placeholder






In addition, I think that the latest Trinity and its changes are responsible for
the contig naming changes.  Its output file, Trinity.fasta, is just the short name 

	comp5_c0_seq2

versus the long name style in the examples 

	comp0_c0_seq1_FPKM_all:1241158.424_FPKM_rel:626919.264_len:482_path:[1032,746,1102]

This difference has led to the modifications of the parsing in 2 of the scripts
and I think that that is causing some problems.

This also effected the write\_results.pl script in another place.  The sequences
weren't included due to this change in name.  The hash then didn't contain any
matching keys.



Also, many of these scripts seem to expect a "standard" sequence name style
of having a trailing /1 or /2 for the respective left and right lanes.
Our data is not like this.  

### TODO 
  Is this really a standard?  Change our sequence names or modify all of the scripts to work with both style?





## RINS modifications ( scripts to be copied in via Makefile )

Modified rins.pl for several reasons.  One major change is to the parameters
passed to Trinity.  It appears that this was written for trinityrnaseq-r20110519.

	Removed ...

		--run_butterfly : no longer needed (or allowed) as is default

		--compatible_path_extension : triggers java exception as is invalid
			Exception in thread "main" java.lang.NullPointerException
				at gnu.getopt.Getopt.checkLongOption(Getopt.java:869)
				at gnu.getopt.Getopt.getopt(Getopt.java:1119)
				at TransAssembly_allProbPaths.main(TransAssembly_allProbPaths.java:193)
			Not valid since trinityrnaseq-r20110519
			( Butterfly/src/src/TransAssembly_allProbPaths.java )

	Changed ...

		--paired_fragment_length changed to --group_pairs_distance
			Although I don't know that that was the right thing to do.

	Added ...

		--JM 1G : or some size
			10G caused a lot of disk activity as I don't have 10G of memory.
			Downsized to 1G and the test_sample run was done in 20 minutes vs 12 hours!


Some file checking as been added.  Many defaults are also now set and are
unnecessary at the command line or in the config file.

Modified blastn\_cleanup.pl and write\_result.pl as the parsing no longer works.
It appears that the field content has changed?


## INSTALLATION

Requirements may include libpng, coreutils, texlive
	sudo port install autoconf automake libpng coreutils texlive


	make 
	make install


This will create ~/RINS\_BASE with all of the scripts and binaries in it.
You will need to modify your PATH to include this.  Feel free to rename it.

setenv PATH ${HOME}/RINS\_BASE/bin:${PATH}


### TODO

Make FASTA copies of our FASTQ data and use it to avoid all the copying and converting.


### TODO

blat can run for quite a long time without ever producing any output.
Make it more verbose. ( using -dots=10000 is helpfulish )


### TODO

Add testing.




### TODO

Merged PathSeq fasta files for all\_bacterial and virus\_all
Blat can't use it as it is too big (>4GB).
Using it to create a blastdb

	makeblastdb -dbtype nucl -in all_bacterial_and_viral.fa -out all_bacterial_and_viral -title all_bacterial_and_viral -parse_seqids

Having trouble. Generated db isn't searchable.  Think I need the parse seqids option...
Yes.  -parse\_seqids is really needed when creating from a fasta file.
Otherwise, new ones are made up which serves no purpose.


### NOTE

Always use -parse\_seqids option when making a blastdb from a fasta file.
Otherwise the sequence names in the fasta file are not used and
new ones are made up.


### TODO

Bowtie2 won't run on my MacPro

	Executing ...
	bowtie2 -p 6 -f /Volumes/cube/working/indexes/hg19 leftlane_not_hg18.fa -S leftlane_not_hg18_hg19.sam --un leftlane_not_hg18_hg19.fa
	bowtie2-align(20745) malloc: *** mmap(size=965771264) failed (error code=12)
	*** error: can't allocate region
	*** set a breakpoint in malloc_error_break to debug
	Out of memory allocating the ebwt[] array for the Bowtie index.  Please try
	again on a computer with more memory.
	Error: Encountered internal Bowtie 2 exception (#1)
	Command: /Users/jakewendt/RINS_BASE/bin/bowtie2-align --wrapper basic-0 -p 6 -f --passthrough /Volumes/cube/working/indexes/hg19 leftlane_not_hg18.fa 

I downloaded the mac binaries and they run fine.  




## SUPER NOTE

Some scripts have the expectation of a trailing lane identifier (either /1 or /2)
Our data is not like this.  It will have something like "\_1:Y:0:CGTACG"
To replace our special lane identifier with the "standard" try something like ...

	cat leftlane.fa | sed 's/_([12]).*$/\/\1/' > relaned.fa 

The underscore is added by the fastq2fasta.pl script so if using fastq try ...

	cat leftlane.fq | sed 's/ ([12]).*$/\/\1/' > relaned.fq

In addition, some scripts, samtools in particular, REQUIRE that the sequence
names DO NOT begin with an @ symbol.

	cat raw.fa | sed 's/^>@/>/' > renamed.fa

These @ symbols are being left behind by RINS' fastq2fasta.pl script

The problem with both of these problems is that they happen quietly in the night.


Trinity's Chrysalis complains about the lane names too.

	util/scaffold_iworm_contigs.pl 
	warning, ignoring read: HWI-ST977:132:C09W8ACXX:7:1101:10003:190695_2:N:0:GTGGCC since cannot decipher if /1 or /2 of a pair.
	Use of uninitialized value $core_acc in string ne at scaffold_iworm_contigs.pl line 64, <$fh> line 1.
	Use of uninitialized value $frag_end in hash element at scaffold_iworm_contigs.pl line 71, <$fh> line 1.




### TODO

Can't make trinity r2012-10-05 on my MacPro or Steve's
It will generate many, seemingly random, ...

	/bin/sh: fork: Resource temporarily unavailable


Solution.  Well, not quite.  These failures appear to occur during the compilation
of the coreutils.  The script actually deals with the failed compilation so
this newer script is not required.  However, Chrysalis has a path hard coded in its 
code so it is expecting a sort command in trinity-plugins/coreutils/bin.
The sole purpose is to allow or deal with a sort command which can do so --parallel.  
sudo port install coreutils will install the same, newer actually, however it
will be gsort.  The contents of /opt/local/libexec/gnubin/ will contain links 
to these "g" utils.  We could install this port, and manually copy gsort to
our path and then this script will notice it and not run the part
that will fail.



### TODO

Change all of the perl shebang lines to "/usr/bin/env perl"







### TODO

Try bwa?

### TODO

Try some other assemblers ....


Price



Velvet
	http://www.ebi.ac.uk/~zerbino/velvet/
	https://github.com/dzerbino/velvet


Ray/RayPlatform
	I compiled and ran this.  Seemed to work, but was taking forever
	so it was killed.  Perhaps give it another go on a big machine.
	http://denovoassembler.sourceforge.net/
	https://github.com/sebhtml/ray.git
	https://github.com/sebhtml/RayPlatform
	(  have "ray" and "RayPlatform" parallel.
	  "ray" has a link to this expected location of "RayPlatform". )

	If it works, add them to our repo as a submodule ...
	git submodule add https://github.com/sebhtml/ray.git
	git submodule add https://github.com/sebhtml/RayPlatform.git

	submodules suck.  They get "dirty" and would need cleaned before commit
	just going to download the code locally.



mira


### Notable Links for Reference

I was researching filtering and repeat masking ...

	http://www.animalgenome.org/bioinfo/resources/manuals/blast2.2.24/
	http://www.girinst.org/server/RepBase/index.php
	http://www.repeatmasker.org/
	http://www.clarkfrancis.com/blast/Blast_Nmer_Masking.html
	http://www.compbio.ox.ac.uk/analysis_tools/BLAST/BLAST_filtering.shtml

Even looked into creating a megablast index

	http://bioinformatics.oxfordjournals.org/content/early/2008/06/23/bioinformatics.btn322.full.pdf




----------
Copyright (c) 2012 [Jake Wendt], released under the MIT license
