# CCLS Sequencing

This is simply a collection of notes and scripts.

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









## New Stuff

I have begun including source code and making modifications to allow
for the complete install of Blat, Blast, Bowtie, Bowtie2, Trinity and RINS
from this repository.  It should be clearly noted that the aforementioned
software packages ARE NOT MINE and belong to those that wrote them.

I created the root Makefile to control the making and installing of these 
packages on my machine and so far this works.

3 blat files required minor modifications.


make clean
make 
make install






## TODO

Add testing.





----------
Copyright (c) 2012 [Jake Wendt], released under the MIT license
