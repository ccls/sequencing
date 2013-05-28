#!/usr/bin/env ruby

require 'bio'
require 'optparse'
require 'ostruct'

o = OpenStruct.new({
	:total_sequences => 0,
	:match_regex     => 'human|sapien'
})

#	nt.fa (2012 version)
#	45032248369 Sep 14  2012 /Volumes/cube/working/indexes/nt.fa
o.total_sequences = 16502268

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB

Usage: #{File.basename($0)} [options] file_name(s)

EOB

	# Define the options, and what they do

	opts.on( '-s', '--total_sequences INTEGER', Integer,
			"Total sequences in input file (used only for reporting)" ) do |s|
		o.total_sequences = s
	end

	opts.on( '-m', '--match STRING', String,
			"Regex to match in sequence name (#{o.match_regex})" ) do |s|
		o.match_regex = s
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
 
##################################################

ARGV.each do |arg|

	puts arg

	if File.exists?(arg)
		inputfile = Bio::FlatFile.auto( arg )
		matched_fasta = File.open("#{arg}.matched",'w')
		matched = 0
		unmatched = 0

		inputfile.each_with_index do |e,i|
			printf "\rProcessing %8d of %8d : matched %8d : unmatched %8d", 
				i, o.total_sequences, matched, unmatched
			if e.definition.match(/#{o.match_regex}/i)
				matched += 1
				matched_fasta.puts e
			else
				unmatched += 1
			end
		end	#	file.each do |sequence|
		matched_fasta.close
		puts
	else	#	if File.exists?(arg)
		puts "#{arg} file not found."
	end	#	if File.exists?(arg)
end


exit;
__END__

inputfile = Bio::FlatFile.auto('nt.fa')
#puts "Counting reads in file"
#> grep "^>" nt.fa | wc -l
# 16502268
#	counting takes WAY WAY too long and why do it each time
#total_sequences = inputfile.count
#puts "Found #{total_sequences}"
#inputfile.rewind

total_sequences = 16502268

human_fasta = File.open('nt_human.fa','w')
human = 0
nonhuman = 0

inputfile.each_with_index do |e,i|
	printf "\rProcessing %8d of %8d : human %8d : non-human %8d", i, total_sequences, human,nonhuman
#	puts e.definition unless e.definition.match(/human|homo\s*sapien|h\.*\s*sapien|home\s*sapien|omo\s*sapien|oma\s*sapien/i)
	if e.definition.match(/human|sapien/i)
		human += 1
		human_fasta.puts e
	else
		nonhuman += 1
	end
end	#	file.each do |sequence|

human_fasta.close

puts

__END__

> bowtie2-build nt_human.fa nt_human
Settings:
  Output files: "nt_human.*.bt2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Max bucket size: default
  Max bucket size, sqrt multiplier: default
  Max bucket size, len divisor: 4
  Difference-cover sample period: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  nt_human.fa
Reading reference sizes
Error: Reference sequence has more than 2^32-1 characters!  Please divide the
reference into batches or chunks of about 3.6 billion characters or less each
and index each independently.
  Time reading reference sizes: 00:01:54
Total time for call to driver() for forward index: 00:01:54
Error: Encountered internal Bowtie 2 exception (#1)
Command: bowtie2-build nt_human.fa nt_human 

Guessing that the "reference sequence" the file and not an individual sequence
6229057055 May 23 15:11 /Volumes/cube/working/indexes/nt_human.fa

Split the file into 3 pieces.  The first was 3685031844 but still worked?



