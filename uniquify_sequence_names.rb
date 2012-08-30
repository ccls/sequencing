#!/usr/bin/env ruby

infilename = "indexes/virus.fa.ORIGINAL"
outfilename = "#{infilename}.out"
sequences = []

File.open(  infilename, 'r' ) { |infile|
File.open( outfilename, 'w' ) { |outfile|
	while( line = infile.gets )
		if line.match(/^>/)
			if sequences.include?(line)
				line.chomp!
				line << "_duplicate\n"	#	if there is a triplicate, will need to run again
			end
			sequences << line
		end
		puts line
		outfile.puts line
	end
} }	#	File.open
