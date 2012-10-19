#!/usr/bin/env ruby
#
# blatoutcandidate.pl loops over the lines in the given psl files and
#   gets the 10ths column ... @HWI-ST281_0133:3:1:6254:2049#0/1
#   removes the trailing "/0" or "/1" leaving @HWI-ST281_0133:3:1:6254:2049#0
#   and adds it to the hash
# Then it loops of the first fasta file 
#   (left hopefully as it writes to blat_out_candidate_leftlane.fa)
#   If the sequence line name is in the hash,
#     prints it and the sequence to the output fasta file
# Then same thing for the second fasta file
#   (right hopefully as it writes to blat_out_candidate_rightlane.fa)
#
# Interesting that both psl files are actually used for each fasta file.
#

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rins_ccls_lib'


#
#	TODO allow for the output filenames option
#


names = {}
(1..(ARGV.length/2)).each do |arg|
	File.open(ARGV[arg-1],'r') do |psl_file|
		5.times{ psl_file.gets }	#	skip header lines (perl version only skips 1???)
		while( line = psl_file.gets )
			data = line.split(/\t/)
			names[data[9].delane_sequence_name] = true
		end
	end
end

if ARGV.length == 4
	#	pair end
	File.open(ARGV[2],'r') { |input|
	File.open('blat_out_candidate_leftlane.fa','w') { |output|
		while( line = input.gets )
			name = line.delane_sequence_name
			if names[name]
				output.puts line
				output.puts input.gets
			else
				input.gets
			end
		end
	} }
	File.open(ARGV[3],'r') { |input|
	File.open('blat_out_candidate_rightlane.fa','w') { |output|
		while( line = input.gets )
			name = line.delane_sequence_name
			if names[name]
				output.puts line
				output.puts input.gets
			else
				input.gets
			end
		end
	} }
elsif ARGV.length == 2
	#	single end
	File.open(ARGV[1],'r') { |input|
	File.open('blat_out_candidate_singlelane.fa','w') { |output|
		while( line = input.gets )
			name = line.delane_sequence_name
			if names[name]
				output.puts line
				output.puts input.gets
			else
				input.gets
			end
		end
	} }
else
	puts "Unexpected number of arguments"
end
