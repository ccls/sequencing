#!/usr/bin/env ruby

usage =<<EOUSAGE

Usage: #{File.basename($0)} non_human_contig.fa non_human_contig_blastn.txt leftlane.psl rightlane.psl results.txt

(Could be 1 or 2 psl files.)

used to clean up data based if the length > 0 and the bit_score/length >= 0.75 and then print out results

Input fasta can now have multi-line sequences.

EOUSAGE

#	show help if arg length is wrong or 
#	used -h or --help (wouldn't be 3 or 6 anyways)
if( ( !ARGV.empty? and ARGV[0].match(/-h/i) ) or
	( ![4,5].include?(ARGV.length) ) )
	puts usage
	exit
end

ARGV[0..-2].each do |f|
	unless File.exists?(f)
		puts "\nFile #{f} not found.\n"
		puts usage
		exit 
	end
end

contig_fa_file     = ARGV[0]			#	first	(non_human_contig.fa)
contig_blastn_file = ARGV[1]			#	second (non_human_contig_blastn.txt)
output_file        = ARGV[-1]			#	last (results.txt)
psl_files          = ARGV[2..-2]	#	third to next to last (leftlane.psl rightlane.psl)
data = {}

File.open(contig_fa_file,'r') do|f|
#	Expecting ....
#	>comp629_c0_seq1 len=320 path=[298:0-60 359:61-166 465:167-207 506:208-235 534:236-259 558:260-271 570:272-319]
#ATAGTTTCTATAGTTTAGTTTATTGCTGTATCATGCTTTCTGGCACGGCAAACTGTCTCCATTGCAAATTTAACAGCTTCTGGGCACCAACTTATTATGACTACTTTCACATAATTACTGTTCTGGCTGCGTTTTCTAGGTTGCCTTGCCAATAAAATGTGCTTCCAAATCTCCACCAAGACACACCTAATCCGGTCGCTGCTTGCTTTCTAGCTTTAATTAATGCAGTTGCTACACGTCTCTTTCTAACTATAATTATAAACTATAATCTAGACAATAATAAATAGGGAGGGACCGAATACGGTGCGACCGAATGGGGT

#	while name_line = f.gets and sequence_line = f.gets
#		#	>comp79_c0_seq1 len=714 path=[692:0-713]
#		name_parts = name_line.chomp.split
#		name   = name_parts[0].gsub(/^>/,'')
#		length = name_parts[1].split(/=/).last
##		puts "#{name}:#{length}"
#		data[name] = {}
#		data[name][:count]     = 0	#	prepare to handle counting
##		data[name][:bitscores] = []	#	prepare to handle multiple bitscores from blastn
#		data[name][:length]    = length
#		data[name][:name_line] = name_line.chomp
#		data[name][:sequence]  = sequence_line.chomp
#	end

	name = ''
	while line = f.gets
		if line =~ /^>/
			#	>comp79_c0_seq1 len=714 path=[692:0-713]
			name_parts = line.chomp.split
			name   = name_parts[0].gsub(/^>/,'')
			length = name_parts[1].split(/=/).last
			#		puts "#{name}:#{length}"
			data[name] = {}
			data[name][:count]     = 0	#	prepare to handle counting
			#		data[name][:bitscores] = []	#	prepare to handle multiple bitscores from blastn
			data[name][:sequence]  = ''	#	prepare to handle multi-line sequence
			data[name][:length]    = length
			data[name][:name_line] = line.chomp
		else
			data[name][:sequence] << line.chomp	#	chomp it as outputing to text file
		end
	end

end

psl_files.each do |psl_file|
	File.open(psl_file,'r') do |f|

##	why 5? (the first 5 lines are the header)
#  1 psLayout version 3
#  2 
#  3 match mis-  rep.  N's Q gap Q gap T gap T gap strand  Q         Q     Q     Q   T         T         T     T   block blockSizes  qStarts  tStarts
#  4       match match     count bases count bases         name      size  start end name      s    ize  start end count
#  5 -------------------------------------------------------------------------------------------    --------------------------------------------------------------------
#  6 100 0 0 0 0 0 0 0 + @HWI-ST281_0133:3:1:11853:3063#0/1  100 0 100 comp621_c0_seq1 1411  18      118 1 100,  0,  18,

		5.times{ f.gets }

		while( line = f.gets )
			line.chomp!
			parts = line.split
			data[parts[13]][:count] += 1
#    my @data = split /\t/, $line;
#    $count{$data[13]} ++; 

		end
	end
end

flag = {}
File.open(output_file, 'w') do |out|
	out.puts "contig_name\tnumber_of_raw_reads_fall_on_this_contig\tnon_human_species\tE-value\tbit_score\tcontig_sequence\n"

	File.open(contig_blastn_file,'r') do |input|
		while line = input.gets
			line.chomp!
	
	#	blastn -outfmt 6 ...
	#	'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
	#	comp9_c0_seq1	gi_109255272_ref_NC_008168.1__Choristoneura_occidentalis_granulovirus,_complete_genome	0.00	0	0	0	0	0	0	0	3e-19	99.0
	
			parts = line.split
			name = parts[0]
			bit_score = parts[-1]
	
			if( ( data[name][:length].to_f > 0 ) && 
					( bit_score.to_f/data[name][:length].to_f >= 0.75) )
	
				pline = "#{parts[0]}\t#{data[parts[0]][:count]}\t#{parts[1]}\t#{parts[10]}\t#{parts[11]}\t#{data[parts[0]][:sequence]}";
	
				#	only print it once
				if( flag[pline].nil? )
					out.puts pline
					flag[pline] = 'not nil'
				end
			end
		end
	end
end
