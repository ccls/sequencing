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

ARGV.each do |filename|
	inputfile     = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)
	sequences = {}
	lanes = []

	puts "Counting sequences"
#	total_sequences = inputfile.count
	total_sequences = `grep '>' #{filename} | wc -l`.to_i
#	will return string with leading spaces and trailing carriage return.  
#=> " 8899236\n"
#	to_i strips it all off and returns just the integer
#=> 8899236

	puts "Found #{total_sequences} sequences"

#	this is probably faster
#grep '>' ~/github_repo/ccls/sequencing/trinity_input_single.fasta  | wc -l
# 8899236
#	yep.  about 10-15x faster

	#	using +1 as log of 2-10 is 1, 11-100 is 2
	max_digits = Math.log10( total_sequences + 1 ).ceil

	inputfile.each_with_index {|e,i|
		printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", i+1, total_sequences
		lane = e.definition.lane.to_i
		lanes.push(lane) unless lanes.include?(lane)
		sequences[e.definition.delane] ||= {}
		sequences[e.definition.delane][lane] = e
	}
	puts	#	mostly for newline after status

	unpaired_sequences = sequences.select{|k,v| v.length == 1 }
	paired_sequences   = sequences.select{|k,v| v.length == 2 }
	confused_sequences = sequences.select{|k,v| v.length > 2 || v.length < 1 }

	puts "Found #{unpaired_sequences.length} unpaired sequences"
	puts "Found #{paired_sequences.length} paired sequences"

	streams = {}
	lanes.each do |lane|
		streams[lane] = File.open("#{root_filename}_#{lane}#{root_extname}",'w')
	end

	total_paired_sequences = paired_sequences.length
	max_digits = Math.log10( total_paired_sequences + 1 ).ceil
	paired_sequences.each_with_index do |kv,i|
		printf "\rWriting sequence %#{max_digits}d of %#{max_digits}d", i+1, total_paired_sequences
		kv[1].each do |lane,sequence|
			streams[lane].puts sequence
		end
	end
	puts	#	for newline after status line

	streams.each{|k,v|v.close}
end

__END__


This version is RIDICULOUSLY slow for big files.

ARGV.each do |filename|
	inputfile     = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)

	sequence_names = Hash.new(0)
	inputfile.each do |entry|
		#	HS1:209:c0v70acxx:6:1101:1633:2228_1:N:0:GATCAG
		#	=> HS1:209:c0v70acxx:6:1101:1633:2228
		#	HS1:209:c0v70acxx:6:1101:1633:2228 1:N:0:GATCAG
		#	=> HS1:209:c0v70acxx:6:1101:1633:2228
		#	HS2:360:D1NTUACXX:8:1101:12120:2042/1
		#	=> HS2:360:D1NTUACXX:8:1101:12120:2042
		#	comp2_c0_seq1 len=295 path=[296:0-50 101:51-172 347:173-294]
		#	=> comp2_c0_seq1 len=295 path=[296:0-50 101:51-172 347:173-294]
		#	
		sequence_names[entry.definition.delane] += 1
	end

	#	only those that were there twice (left and right)
	sequence_names.delete_if{|name,count| count < 2}

	#	convert to array of just the names
	sequence_names = sequence_names.collect{|name,count| name }

	#	open the output files and keep in order in this array
	streams = [File.open("#{root_filename}_1#{root_extname}",'w'),
		File.open("#{root_filename}_2#{root_extname}",'w')]

#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#	streams = [File.open("#{root_filename}_0#{root_extname}",'w'),
#		File.open("#{root_filename}_1#{root_extname}",'w'),
#		File.open("#{root_filename}_2#{root_extname}",'w')]

	inputfile.rewind	#	already read the file, so need to rewind
	inputfile.each do |entry|
		if sequence_names.include?(entry.definition.delane)
			#	if matches, write to the correct file
			streams[entry.definition.lane.to_i - 1].puts entry
#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#			streams[entry.definition.lane.to_i].puts entry
		end
	end
	streams.each{|f|f.close}

end	#	ARGV.each do |filename|
__END__
