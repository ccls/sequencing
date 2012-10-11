#!/usr/bin/env ruby

#
#	This is my version of the cleanup script.
#	As the trinity output has changed, we need to get
#	the contig length from the fasta file as it isn't
#	in the blastn output.
#
#	AHHH!  This is why the modify_trinity_output is called
#	It puts the sequence on one line. 

# blastn_cleanup.rb human_contig.txt Trinity.fasta clean_blastn.fa 0.8


#	read the blastn output
#	read the fasta
#	THEN
#	loop through the combined hash data

#	data = {}
#	data['comp430_c0_seq1'] = {
#		:length(from_fasta)
#		:bitscore(from blastn)	actually bitscores or max_bitscore
#		:sequence
#	}


raise "BAD DOG.  Expected 4 args." if ARGV.length != 4

input = ARGV[0]            #	human_contig.txt
fasta_input = ARGV[1]      #	Trinity.fasta
output = ARGV[2]           #	clean_blastn.fa
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
		if( ((max_bitscore.to_i)/(v[:length].to_i)) >= similarity_thrd.to_i )
			f.puts v[:name_line]
			f.puts v[:seq_line]
		end
	end
end
