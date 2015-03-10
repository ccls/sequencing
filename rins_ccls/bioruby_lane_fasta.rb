#!/usr/bin/env ruby

require 'bio'	#	bioruby gem
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'

#	read a fasta file
#	find the reads that have both a left AND right lane (1 and 2)
#		(sometimes this is 0 and 1)
#	

if ARGV.empty?
	puts "#{$0} filename(s)\n\nAt least one input fasta file name is required.\n\n"
	exit
end

#	Runs out of memory on n0 with 68GB fasta file

ARGV.each do |filename|
	inputfile     = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)
	lanes = [1,2]	#[]

	#command = "grep '^>' #{filename} | wc -l"
	#	grep has a built-in counter which is faster
	command = "grep -c '^>' #{filename}"
	puts "Counting sequences with ..."
	puts command

	#	`` will return string with leading spaces and trailing carriage return.  
	#=> " 8899236\n"
	#	to_i strips it all off and returns just the integer. Awesome.
	##=> 8899236
	total_sequences = `#{command}`.to_i

	puts "Found #{total_sequences} sequences"

	#	using +1 as log of 2-10 is 1, 11-100 is 2
	max_digits = Math.log10( total_sequences + 1 ).ceil


	#	This is blunt and effective, I think, but it all happens behind the scenes
	#	It would be nice to make this more verbose.
	#	It is substantially easier on memory though as it doesn't load
	#	the entire file into memory.
	puts "Scanning #{filename} for paired sequences (this can take a while) with ..."


#	this expects sequence names to end with /1 or /2
	command = "grep '^>' #{filename} | awk -F/ '{print $1}' | sort | uniq -d | sed 's/^>//'"
#	this expects a space then lane-related and other irrelevant stuff
#	command = "grep '^>' #{filename} | awk '{print $1}' | sort | uniq -d | sed 's/^>//'"
#	As the purpose of this is to prep files for use with Trinity and
#	Trinity requires the /1 or /2 laning, need to use the /1 or /2 laning.


	puts command
	paired_sequences_a = `#{command}`.chomp.split
	paired_sequences = {}
	paired_sequences_a.each{|i| paired_sequences[i] = 1 }
	paired_sequences_a = nil

	puts "Found #{paired_sequences.length} paired sequences"

	puts "Opening output files."
	streams = {}
	lanes.each do |lane|
		streams[lane] = File.open("#{root_filename}_#{lane}#{root_extname}",'w')
	end

	puts "Writing to output files."
	total_paired_sequences = paired_sequences.length



	index=0
	inputfile.each do |entry|
#	inputfile.each_with_index do |entry,index|

		printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", index+1, total_sequences
		#	WAY WAY WAY TOO SLOW.  OMG
		#		if paired_sequences.include?(entry.definition.delane)
		#	Hash#has_key? IS WAY SUPER FASTER COMPARED TO Array#includes?
		#	I would've thunk it the opposite.
		if paired_sequences.has_key?(entry.definition.delane)
			#
			#	Add some verbosity here?
			#

			#	if matches, write to the correct file
			streams[entry.definition.lane.to_i].puts entry
#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#			streams[entry.definition.lane.to_i].puts entry
		end


		index+=1	#	manually increment
		nil	#	reduce memory usage
	end
	puts
	puts "Closing output files."
	streams.each{|k,v|v.close}
end

__END__

