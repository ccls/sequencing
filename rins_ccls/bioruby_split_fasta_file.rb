#!/usr/bin/env ruby

require 'bio'	#	bioruby gem
require 'optparse'
require 'ostruct'

puts "Running #{$0}"
o = OpenStruct.new({
	:verbose => false
})

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB

Usage: #{File.basename($0)} [options] file_name(s)

EOB

	# Define the options, and what they do
	opts.on( '-v', '--verbose', 
			"Print out more info" ) do |s|
		puts "Verbosity is on."
		o.verbose = true
	end


	opts.on( '-r', '--read_count INTEGER', Integer,
			"Max number of reads per output file" ) do |s|
		o.read_count = s
	end

	opts.on( '-f', '--file_count INTEGER', Integer,
			"Max number of files per output" ) do |s|
		o.file_count = s
	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end

	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
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

if ARGV.empty?
	puts optparse	#	Basically display the command line help
	puts
	puts "At least one file name is required."
	puts
	exit
end

if o.respond_to?(:read_count) and o.respond_to?(:file_count)
	puts optparse	#	Basically display the command line help
	puts
	puts "Cannot specify both read_count and file_count."
	puts
	exit
end

#	if neither read_count nor file_count are given
#	use a default of 10 reads per output file?
unless o.respond_to?(:read_count) or o.respond_to?(:file_count)
	puts
	puts "No count given.  Using read_count of 10."
	puts
	o.read_count = 10
end

ARGV.each do |filename|

	inputfile = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)

	#	Bio::FastaFormat or Bio::Fastq (don't know why not FastqFormat)
	puts "Detected #{inputfile.dbclass} format" if o.verbose	

	total_sequences = inputfile.count
	inputfile.rewind

	puts "Counted #{total_sequences} sequences in file." if o.verbose

	#	MUST use a to_f for integer division or will get integer out (floor not ceil)

	num_files = o.file_count || ( total_sequences.to_f / o.read_count ).ceil
	#	needed so know how many leading zeroes to use
	puts "Number of output files to create ... #{num_files}" if o.verbose

	num_reads_per_file = o.read_count || ( total_sequences.to_f / o.file_count ).ceil
	puts "Max number of reads per output file ... #{num_reads_per_file}" if o.verbose

	#	TODO how about an outdir?

	#	as there is no each_slice_with_index
	slice_counter = 0
	inputfile.each_slice(num_reads_per_file) do |slice|
		slice_counter += 1
		outfilename = "#{root_filename}_#{sprintf(
				"%0#{Math.log10(num_files).ceil}d",slice_counter)}#{root_extname}"
		puts "Writing to #{outfilename}"
		File.open(outfilename,'w'){|f| f.puts slice }
	end	#	file.each do |sequence|

end	#	ARGV.each do |filename|

__END__
