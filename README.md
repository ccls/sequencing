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



## Trinity modifications ( scripts to be copied in after via Makefile )

SAM\_filter\_out\_unmapped\_reads.pl was modified due to some test
data triggering a divide by zero error.


## RINS modifications ( scripts to be copied in after via Makefile )

Modified rins.pl for several reasons.  One major change is to the parameters
passed to Trinity.  It appears that this was written for trinityrnaseq-r20110519.

Some file checking as been added.  Many defaults are also now set and are
unnecessary at the command line or in the config file.

Modified blastn\_cleanup.pl and write\_result.pl as the parsing no longer works.
It appears that the field content has changed?


## INSTALLATION

make clean
make 
make install


This will create ~/RINS\_BIN with all of the scripts and binaries in it.
You will need to modify your PATH to include this.  Feel free to rename it.

setenv PATH ${HOME}/RINS\_BIN:${PATH}


## TODO

Add testing.


----------
Copyright (c) 2012 [Jake Wendt], released under the MIT license
