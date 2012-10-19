#!/usr/bin/env ruby

require 'optparse'

# This hash will hold all of the options
# parsed from the command-line by
# OptionParser.
options = {
	:verbose => false
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	opts.banner = "\nUsage: #{File.basename($0)} [options] -i in_file -o out_file\n\n"

	#	Define the options, and what they do

	opts.on( '-i', '--in FILE', 'tab delimited txt infile (required)' ) do |file|
		options[:in] = file
	end

	opts.on( '-o', '--out FILE', 'tab delimited txt outfile (required)' ) do |file|
		options[:out] = file
	end

	opts.on( '-v', '--verbose', 'Output more information' ) do
		options[:verbose] = true
	end

	# This displays the help screen, all programs are assumed to have this option.
	#	Add extra "\n" to last option for aesthetics.
	opts.on( '-h', '--help', 'Display this help screen',"\n") do
		puts opts
		exit
	end
end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options. What's left is the list of files to resize.
optparse.parse!

#	Infile and outfile required
unless options[:in] and options[:out]
	puts optparse	#	Basically display the command line help
	exit
end

unless File.exists?(options[:in])
	puts "File #{options[:in]} not found." 
	exit
end
 
File.open(  options[:in], 'r' ) { |infile|
File.open( options[:out], 'w' ) { |outfile|

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
