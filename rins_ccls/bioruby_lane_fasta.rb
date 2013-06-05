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
