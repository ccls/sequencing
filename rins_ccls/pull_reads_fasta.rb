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

usage =<<EOUSAGE

Usage: #{File.basename($0)} name_files fasta_files output_files

The fasta input files are currently expected to have the entire sequence on a single line.

pull_reads_fasta.rb joins the sequence lists from the name_files, then 'de-lanes' the list, removing the left/right, 0/1 or 1/2 suffixes, and uniquing them in order to keep both lanes in sync, then selects those sequences from the inputs and placing them in the outputs.

EOUSAGE

#	show help if arg length is wrong or 
#	used -h or --help (wouldn't be 3 or 6 anyways)
if( ( !ARGV.empty? and ARGV[0].match(/-h/i) ) or
	( ![3,6].include?(ARGV.length) ) )
	puts usage
	exit
end

#	Check all input files ... 0,1 or 0,1,2,3
#ARGV[0..(ARGV.length-(ARGV.length/3)-1)].each do |f|
#	or ...
ARGV[0..(-(ARGV.length/3)-1)].each do |f|
	unless File.exists?(f) 
		puts "\nFile #{f} not found.\n"
		puts usage
		exit 
	end
end

#
#	This is very similar to blatoutcandidate.rb
#

def merge_names(*args)
	names = {}
	args.each do |arg|
		File.open(arg,'r') do |f|
			while line = f.gets do
				names[line.delane_sequence_name] = true
			end
		end
	end
	names
end

if ARGV.length == 6
	#	pair_end
	names = merge_names(ARGV[0],ARGV[1])
	pull_reads(names,ARGV[2],ARGV[4])
	pull_reads(names,ARGV[3],ARGV[5])
elsif ARGV.length == 3
	#	single_end
	names = merge_names(ARGV[0])
	pull_reads(names,ARGV[1],ARGV[2])
else	#	shouldn't happen with above check
	puts
	puts "Unexpected number of arguments ( 3 or 6 )"
	puts
end
