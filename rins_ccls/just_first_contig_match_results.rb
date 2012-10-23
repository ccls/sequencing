#!/usr/bin/env ruby

#	Infile and outfile required
unless ARGV.length == 2
	puts "\nUsage: #{File.basename($0)} in_file out_file\n\n" <<
		"Outputs just the first reference to each contig name.\n" <<
		"The in file is expected to be LIKE blastn outfmt 6.\n\n"
	exit
end

unless File.exists?(ARGV[0])
	puts "File #{ARGV[0]} not found." 
	exit
end
 
File.open( ARGV[0], 'r' ) { |infile|
File.open( ARGV[1], 'w' ) { |outfile|

#	first line better be a header line
#	   contig_name^Inumber_of_raw_reads_fall_on_this_contig^Inon_human_species^IE-value^Ibit_scor     e^Icontig_sequence$

	line = infile.gets
	outfile.puts line

	contig_names = {}

	while( line = infile.gets )
#	  comp155_c0_seq1^I9^Igi|115286659|gb|CP000443.1|^I3e-149^I 536^IGCGGCGTGTTCTTCCTCGGCTTCGGTG     GTTCGCGTTGGACTTCCGACATGGCGATCAACTACTACAAGGTCGTCCTGGGCGTGGCCGCGCAGCTCTTCGCAATGGTGCTCCTGGTGG     GCATCGGCAAGACCTTCCTCGATGACTACTACGCGCGCATGAGCGCCGGCATCAGCCTCAAGGAAATGGGTGTGATGCTGATCGTCGTCA     TCATCCTTCTGGCGTTGACCAACAAGATTCCGCCGCTCATCGCCGGGATCATCACCGGCGCGAGCGTGGGCGGTGCCGGCATCGGTCAGT     TCGGTGCAGGCGCTGCGCT$

		parts = line.split("\t")
		outfile.puts line unless contig_names[parts[0]]
		contig_names[parts[0]] = true

	end
} }	#	File.open
