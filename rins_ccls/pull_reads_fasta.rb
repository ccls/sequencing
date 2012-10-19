#!/usr/bin/env ruby
#
#	The simplified and rubified version of the rins perl script
#
#	There is no standard sequence naming convention so the perl script
#	would not correctly parse our fallon or targgetted data.
#	Rather than modify the perl script, create the ruby script.
#
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'rins_ccls_lib'

def merge_names(*args)
	names = {}
	args.each do |arg|
		File.open(arg,'r') do |f|
			while line = f.gets do
				names[line.chomp] = true
			end
		end
	end
	names
end

def pull_reads( names, input_fasta, output_fasta)
	File.open( input_fasta,'r') { |input|
	File.open(output_fasta,'w') { |output|
		while line = input.gets do
    	if( names[line.delane_sequence_name] )
      	output.puts line
      	output.puts input.gets
    	else
      	input.gets
    	end
		end
	} }
end

if ARGV.length == 6
	#	pair_end
	left_names_in   = ARGV[0]
	right_names_in  = ARGV[1]
	left_fasta_in   = ARGV[2]
	right_fasta_in  = ARGV[3]
	left_fasta_out  = ARGV[4]
	right_fasta_out = ARGV[5]
	names = merge_names(left_names_in,right_names_in)
	pull_reads(names, left_fasta_in, left_fasta_out)
	pull_reads(names,right_fasta_in,right_fasta_out)
elsif ARGV.length == 3
	#	single_end
	single_names_in   = ARGV[0]
	single_fasta_in   = ARGV[1]
	single_fasta_out  = ARGV[2]
	names = merge_names(single_names_in)
	pull_reads(names,single_fasta_in,single_fasta_out)
else
	puts "Unexpected number of arguments"
end
