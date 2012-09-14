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

[Bowtie 0.12.8 and Bowtie 2.0.0 beta 7](http://bowtie-bio.sourceforge.net/)

[Trinity 2012-06-08](http://trinityrnaseq.sourceforge.net) 
Not including the ~50MB of sample data



## Indexes Used

Manually downloaded, but then for some reason couldn't actually use?
ftp://ftp.ncbi.nlm.nih.gov/blast/db/
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz


So used update\_blastdb.pl to get them

	update_blastdb.pl --decompress nt
	update_blastdb.pl --decompress nr

This will downloaded them to whereever you are so make sure you have ~50GB of space!


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



## Trinity modifications ( scripts to be copied in via Makefile )

SAM\_filter\_out\_unmapped\_reads.pl was modified due to some test
data triggering a divide by zero error.

Add an empty target to a Makefile.in.
I doesn't seem needed on my Mac Pro, just my MacBook Pro?

	> trinityrnaseq_r2012-06-08/trinity-plugins/jellyfish-1.1.5/Makefile.in 

		doc/jellyfish.man:
		#      this is just an empty placeholder




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

	make 
	make install


This will create ~/RINS\_BASE with all of the scripts and binaries in it.
You will need to modify your PATH to include this.  Feel free to rename it.

setenv PATH ${PATH}:${HOME}/RINS\_BASE/bin


## TODO

Add testing.


----------
Copyright (c) 2012 [Jake Wendt], released under the MIT license
