#!/usr/bin/env ruby

#
#	This is my version of the cleanup script.
#	The purpose of it is to basically filter out the
#	input fasta based on the bitscores in the contig
#	text file.
#
#	As the trinity output has changed, we need to get
#	the contig length from the fasta file as it isn't
#	in the blastn output.
#

#	read the blastn output
#	read the fasta
#	THEN
#	loop through the combined hash data

#
#	TODO allow multi-line sequences and deprecate remove modify_trinity_output.pl
#

usage =<<EOUSAGE

Usage: #{File.basename($0)} [options] contigs_input fasta_input fasta_output similarity_thrd

The contigs_input file is expected to be in blastn outfmt 6.

The fasta_input file is currently expected to have the entire sequence on a single line.  This is why 'modify_trinity_output.pl' is used.

What exactly is 'similarity_thrd'?

EOUSAGE

#	show help if arg lenght is wrong or used -h or --help (won't be 4 so would work anyways)
if( ( ARGV.length != 4 ) ||
	( !ARGV.empty? and ARGV[0].match(/-h/i) ) )
	puts usage
	exit
end

ARGV[0,1].each do |f|
	unless File.exists?(f)
		puts "\nFile #{f} not found.\n\n" 
		puts usage
		exit
	end
end

input           = ARGV[0]  #	human_contig.txt
fasta_input     = ARGV[1]  #	Trinity.fasta
output          = ARGV[2]  #	clean_blastn.fa
similarity_thrd = ARGV[3]  #	0.8

data = {}

File.open(fasta_input,'r') do|f|
#	Expecting ....
#	>comp629_c0_seq1 len=320 path=[298:0-60 359:61-166 465:167-207 506:208-235 534:236-259 558:260-271 570:272-319]
#ATAGTTTCTATAGTTTAGTTTATTGCTGTATCATGCTTTCTGGCACGGCAAACTGTCTCCATTGCAAATTTAACAGCTTCTGGGCACCAACTTATTATGACTACTTTCACATAATTACTGTTCTGGCTGCGTTTTCTAGGTTGCCTTGCCAATAAAATGTGCTTCCAAATCTCCACCAAGACACACCTAATCCGGTCGCTGCTTGCTTTCTAGCTTTAATTAATGCAGTTGCTACACGTCTCTTTCTAACTATAATTATAAACTATAATCTAGACAATAATAAATAGGGAGGGACCGAATACGGTGCGACCGAATGGGGT
	while name_line = f.gets and sequence_line = f.gets
		#	>comp79_c0_seq1 len=714 path=[692:0-713]
		name_parts = name_line.chomp.split
		name   = name_parts[0].gsub(/^>/,'')
		length = name_parts[1].split(/=/).last
#		puts "#{name}:#{length}"
		data[name] = {}
		data[name][:bitscores] = []	#	prepare to handle multiple bitscores from blastn
		data[name][:length]    = length
		data[name][:name_line] = name_line.chomp
		data[name][:seq_line]  = sequence_line.chomp
	end
end

File.open(input,'r') do|f|
#	comp629_c0_seq1	gi|9627257|ref|NC_001576.1|	100.00	65	0	0	256	320	65	5e-25	 121
# blastn -outfmt 6 ...
# 'qseqid         sseqid pident length mismatch gapopen qstart qend 
# comp430_c0_seq1 chrX   97.24   398      10      1       1    398  
# sstart       send    evalue bitscore'
# 108184401  108184005   0.0    673
	while line = f.gets
		parts = line.chomp.split
		data[parts[0]][:bitscores].push parts[-1]
	end
end

File.open(output,'w') do |f|
	data.each_pair do |k,v|
		max_bitscore = v[:bitscores].max
		puts "#{k}..."
		#	VERY POSSIBLY EMPTY (no blastn matches)
		#puts "  bitscores    : #{v[:bitscores].join(',')}"	#	sometimes quite a few

		puts "  max bitscore : #{max_bitscore}"	#	VERY POSSIBLY BLANK
		#	if there are no bitscores, then max is nil
		#	if there are no bitscores, does that mean that
		puts "  length       : #{v[:length]}"

		percent = ((max_bitscore.to_f)/(v[:length].to_f))

		#	low value means less likely
		puts "  percent      : #{percent}"
		puts "  similarity   : #{similarity_thrd}"
		
		keep = percent <  similarity_thrd.to_f
		puts "  ... #{(keep)? 'keeping' : 'cleaning'}"

		#	original script did a (flag if >= and then print if !flagged)
		#	because compared multiple times. Since using max and only
		#	comparing once can just use < and not the !
		#		if(! ((max_bitscore.to_f)/(v[:length].to_f)) >= similarity_thrd.to_f )

		if( keep )
			f.puts v[:name_line]
			f.puts v[:seq_line]
		end
	end
end
