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
require 'sequencing_lib'

#	blatoutcandidate.rb psl_files fasta_files output_files

usage =<<EOUSAGE

Usage: #{File.basename($0)} psl_files fasta_files output_files

#{File.basename($0)} joins the sequence lists from the psl_files, then 'de-lanes' the list, removing the left/right, 0/1 or 1/2 suffixes, and uniquing them in order to keep both lanes in sync, then selects those sequences from the inputs and placing then in the outputs.

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

names = {}
(1..(ARGV.length/3)).each do |arg|
	File.open(ARGV[arg-1],'r') do |psl_file|
		5.times{ psl_file.gets }	#	skip header lines (perl version only skips 1???)
		while( line = psl_file.gets )
			data = line.split(/\t/)
			names[data[9].delane_sequence_name] = true
		end
	end
end

if ARGV.length == 6
	#	pair end
	pull_reads(names,ARGV[2],ARGV[4])
	pull_reads(names,ARGV[3],ARGV[5])
elsif ARGV.length == 3
	#	single end
	pull_reads(names,ARGV[1],ARGV[2])
else	#	shouldn't happen with above checks
	puts
	puts "Unexpected number of arguments ( 3 or 6 )"
	puts
end
