#!/usr/bin/env ruby
#
#	The simplified and rubified version of the rins perl script
#
#	There is no standard sequence naming convention so the perl script
#	would not correctly parse our fallon or targgetted data.
#	Rather than modify the perl script, create the ruby script.
#
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'

puts
puts "pull_reads_fasta.rb name_files fasta_files output_files"
puts
puts "pull_reads_fasta.rb joins the sequence lists from the name_files,"
puts "\tthen 'de-lanes' the list, removing the left/right, 0/1 or 1/2 suffixes,"
puts "\tand uniquing them in order to keep both lanes in sync,"
puts "\tthen selects those sequences from the inputs and placing then in the outputs."
puts

#
#	This is very similar to blatoutcandidate.rb
#

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
	puts
	puts "Unexpected number of arguments ( 3 or 6 )"
	puts
end
