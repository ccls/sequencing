#!/usr/bin/env ruby

require 'optparse'

# This hash will hold all of the options
# parsed from the command-line by
# OptionParser.
options = {
	:suffix  => 'duplicate',
	:verbose => false
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	opts.banner = "Usage: #{$0} [options] fa_in_file1 fa_in_file2 ...\n" <<
		"Loop over all sequence names and append '#{options[:suffix]}_#' to duplicates."

	# Define the options, and what they do

	opts.on( '-s', '--suffix STRING', 
			"duplicate sequence suffix (default '#{options[:suffix]}_#')" ) do |s|
		options[:suffix] = s
	end

	opts.on( '-v', '--verbose', 'Output more information' ) do
		options[:verbose] = true
	end

	# This displays the help screen, all programs are assumed to have this option.
	opts.on( '-h', '--help', 'Display this screen' ) do
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
 
ARGV.each do |infilename|
	if File.exists? infilename
		puts "Processing #{infilename}"
	else
		puts "#{infilename} not found. Skipping."
		next
	end
#	infilename = "indexes/virus.fa.ORIGINAL"
	outfilename = "#{infilename}.#{$0}.out"
	sequences = []

	File.open(  infilename, 'r' ) { |infile|
	File.open( outfilename, 'w' ) { |outfile|
		while( line = infile.gets )
			if line.match(/^>/)
				puts "Found sequence name."
				puts line
				i=0
				while sequences.include?(line)
					line.chomp!
					line << "_#{options[:suffix]}_#{i}\n"
					puts "Duplicate sequence name. Renaming."
					puts line
					i+=1
				end
				sequences << line
			end
			puts line if options[:verbose]
			outfile.puts line
		end
	} }	#	File.open
end

