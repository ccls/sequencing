#!/usr/bin/env ruby

#	It is sad that the original script is called sam2names, but it is not true.
#	The script only prints the names that have a third column of "*"

#@HWI-ST281_0133:3:1:10412:2250#0/2      0       chr4    83513690        255     100M    * 0       0       CGCCTCCAGACGCGGTTCCGCCCCCGGTCCCGGCTCCGCTTCCCGCCGCCGCCGCTGCCCCCTGTGTCGCCGCCACCATGGCTCCCCCCTGCTCGCCCGC    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    XA:i:1  MD:Z:0T85T13    NM:i:2
#@HWI-ST281_0133:3:1:13342:2066#0/2      4       *       0       0       *       *       0 0       CGGCATTCCAGCGGAACCGCTCGTCCGAGCCTCGGCATTCCTGGGGGCCCCCCCCTTTAAAAAAAAAAAAAAAAAGACGGGGAAGGGAAAGAGGTGAGCA    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    XM:i:0

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'

usage =<<EOUSAGE

Usage: #{File.basename($0)} sam_in_file names_out_file

It is sad that the script is called sam2names, because it is not accurate.

The script only prints the names that have a third column of "*"

EOUSAGE

#	show help if arg length is wrong or 
#	used -h or --help (wouldn't be 3 or 6 anyways)
if( ( !ARGV.empty? and ARGV[0].match(/-h/i) ) or
	( ARGV.length != 2 ) )
	puts usage
	exit
end

unless File.exists?(ARGV[0]) 
	puts "\nFile #{ARGV[0]} not found.\n"
	puts usage
	exit 
end

#
#	This is very similar to blatoutcandidate.rb
File.open( ARGV[0], 'r' ) { |input|   #	sam file
File.open( ARGV[1], 'w' ) { |output|  #	names file
	while( line = input.gets )
		data = line.chomp.split(/\t/)
		if (data[2] == "*")
			output.puts data[0].delane_sequence_name
		end
	end
} }
