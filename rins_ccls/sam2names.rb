#!/usr/bin/env ruby

#	It is sad that the original script is called sam2names, but it is not true.
#	The script only prints the names that have a third column of "*"

#@HWI-ST281_0133:3:1:10412:2250#0/2      0       chr4    83513690        255     100M    * 0       0       CGCCTCCAGACGCGGTTCCGCCCCCGGTCCCGGCTCCGCTTCCCGCCGCCGCCGCTGCCCCCTGTGTCGCCGCCACCATGGCTCCCCCCTGCTCGCCCGC    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    XA:i:1  MD:Z:0T85T13    NM:i:2
#@HWI-ST281_0133:3:1:13342:2066#0/2      4       *       0       0       *       *       0 0       CGGCATTCCAGCGGAACCGCTCGTCCGAGCCTCGGCATTCCTGGGGGCCCCCCCCTTTAAAAAAAAAAAAAAAAAGACGGGGAAGGGAAAGAGGTGAGCA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    XM:i:0

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rins_ccls_lib'

samfile  = ARGV[0];
namefile = ARGV[1];

File.open( samfile,'r') { |input|
File.open(namefile,'w') { |output|
	while( line = input.gets )
		line.chomp!
		data = line.split(/\t/)
		if (data[2] == "*")
			output.puts data[0].delane_sequence_name
		end
	end
} }
