#!/usr/bin/env ruby

require 'bio'	#	bioruby gem
require 'optparse'

number_of_random_reads = 1000

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB

Usage: #{File.basename($0)} [options] input_filename

EOB

	# Define the options, and what they do

	opts.on( '-n', '--number_of_random_reads INTEGER', Integer,
			"Number of random reads to create (default #{number_of_random_reads})" ) do |int|
		number_of_random_reads = int
	end

#	opts.on( '-h', 'Display simple help screen' ) do
#		puts opts
#		puts
#		exit
#	end
#
#	opts.on( '--help', 'Display thorough help screen' ) do
#		puts opts
#		puts
##puts <<EOB
##
##Usage: #{File.basename($0)} --output_suffix "darkdev.fastq.bowtie2.trinity20120608"
##
##Loops through all of the bowtie_human_indexes in the config.yml
##stripping out all of the matches using bowtie2's --un-conc option.
##
##Trinity attempts to assemble the resultant file.
##
##Those results are then blastn'd.
##
##EOB
#		exit
#	end
end

# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options.
optparse.parse!

if ARGV.empty?
	puts optparse   #       Basically display the command line help
	puts
	puts "Input file name is required."
	puts
	exit
end


filename  = ARGV[0]
#	could use a while and loop through all given files




inputfile     = Bio::FlatFile.auto(filename)
root_extname  = File.extname(filename)
root_filename = File.basename(filename,root_extname)

command = "grep '^>' #{filename} | wc -l"
puts "Counting sequences with ..."
puts command
total_sequences = `#{command}`.to_i
#		#	will return string with leading spaces and trailing carriage return.  
#		#=> " 8899236\n"
#		#	to_i strips it all off and returns just the integer. Awesome.
#		##=> 8899236
#
puts "Found #{total_sequences} sequences"
#	using +1 as log of 2-10 is 1, 11-100 is 2
max_digits = Math.log10( total_sequences + 1 ).ceil


puts "\nSelecting #{number_of_random_reads} reads"

if( number_of_random_reads > total_sequences )
	puts
	puts "More reads requested than exist."
	puts
	exit
end


random_read_indexes=(0..total_sequences-1).to_a.sample(number_of_random_reads)

puts "Opening output file."
outputfile = File.open("#{root_filename}.select.#{root_extname}",'w')

inputfile.each_with_index do |entry,index|
	printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", index+1, total_sequences
	if random_read_indexes.include?(index)
		outputfile.puts entry
	end
end
puts
puts "Closing output files."
outputfile.close




#	end while

__END__
