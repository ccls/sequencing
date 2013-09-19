#!/usr/bin/env ruby

#require 'bio'
require 'optparse'
#require 'ostruct'

number_of_random_reads = 1000
length_of_random_reads = 100

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB

Usage: #{File.basename($0)} [options] [optional output filename]

EOB

	# Define the options, and what they do

	opts.on( '-n', '--number_of_random_reads INTEGER', Integer,
			"Number of random reads to create (default #{number_of_random_reads})" ) do |int|
		number_of_random_reads = int
	end

	opts.on( '-l', '--length_of_random_reads INTEGER', Integer,
			"Length of random reads created (default #{length_of_random_reads})" ) do |int|
		length_of_random_reads = int
	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end

	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
		puts
#puts <<EOB
#
#Usage: #{File.basename($0)} --output_suffix "darkdev.fastq.bowtie2.trinity20120608"
#
#Loops through all of the bowtie_human_indexes in the config.yml
#stripping out all of the matches using bowtie2's --un-conc option.
#
#Trinity attempts to assemble the resultant file.
#
#Those results are then blastn'd.
#
#EOB
		exit
	end
end

# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options.
optparse.parse!
 



output = ( ARGV.empty? ) ? STDOUT : File.open(ARGV.shift,'w')



chars = ['A','C','G','T']


digits = Math.log10(number_of_random_reads).ceil

number_of_random_reads.times.each do |i|

	output.printf ">random_%0#{digits}d\n",i
	output.puts length_of_random_reads.times.collect{|i| chars[rand(chars.length)] }.join()

end

output.close unless output == STDOUT
