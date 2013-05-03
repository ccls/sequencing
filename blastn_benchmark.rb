#!/usr/bin/env ruby

require 'csv'

time = Time.now.strftime("%Y%m%d%H%M%S")

c = CSV.open("blastn_times_#{time}.csv",'w')
c << %w( iteration fasta_file evalue num_threads strand num_descriptions num_alignments max_target_seqs word_size ungapped outfmt real_time user_time sys_time )

#	       82.36 real        74.28 user         7.33 sys

( 1..10 ).each do |iteration|

#		blast_input_1_read_201b.fasta
#		blast_input_1_read_201c.fasta

#%w( blast_input_1_read_201a.fasta
#		blast_input_3_reads_201.fasta 
%w(	blast_input_1_read_854.fasta
		blast_input_1_read_895.fasta
		blast_input_1_read_957.fasta
		blast_input_1_read_1754.fasta
		blast_input_16_reads.fasta
			).each do |fasta|
#%w( blast_input_1_read_201a.fasta
#			).each do |fasta|

#(1..4).each do |num_threads|
(1..2).each do |num_threads|
%w( both plus minus ).each do |strand|
[10].each do |num_descriptions| #[10,100,0,1].each do |num_descriptions|
[10].each do |num_alignments| #[10,100,0,1].each do |num_alignments|
[10].each do |max_target_seqs| #[10,100,1].each do |max_target_seqs|
(8..12).each do |word_size| #[4,8,12].each do |word_size|
[false,true].each do |ungapped|
[0].each do |outfmt| #[0,5].each do |outfmt|
%w( 1e-10 1e-50 ).each do |evalue| #%w( 1e-1 1e-10 1e-50 1e-100 ).each do |evalue|

	
#	2.2.27 always raises error with num_threads other than 1 on my macpro?
#	I think that I need to set some flag for multithreading at compilation.
#Error: NCBI C++ Exception:
#    "/Users/jakewendt/github_repo/ccls/sequencing/ncbi-blast-2.2.27+-src/c++/src/corelib/ncbiobj.cpp", line 689: Critical: ncbi::CObject::ThrowNullPointerException() - Attempt to access NULL pointer.
#
#	2.2.25 seems ok


	command =  "/usr/bin/time blastn -db /Users/jakewendt/blast_indexes/nt "
	command << "-num_threads #{num_threads} "
	command << "-query #{fasta} "
	command << "-evalue #{evalue} "	#1e-40 "	#	doesn't seem to make a difference

	command << "-use_index false "	##{use_index}"	#	true or false
	command << "-strand #{strand} "	#	both, plus, minus
	command << "-num_descriptions #{num_descriptions} "	#		>= 0	(default 500)
	command << "-num_alignments #{num_alignments} "	#	>= 0	(default 250)
	command << "-max_target_seqs #{max_target_seqs} "	#	>= 1 (default is ????)
	command << "-word_size #{word_size} "	#		>= 4 (default is ????)	# seems to cause blastn to hang
	command << "-ungapped " if ungapped		#	no value so will need a condition test 
	command << "-outfmt #{outfmt} "	#		default 0, xml = 5


	command << "-out /dev/null "
	command << "2>&1"						#	need this to get the output from time

#		command = "/usr/bin/time echo 'hello world' > /dev/null 2>&1"
#		command = "/usr/bin/time sleep 1 2>&1 > /dev/null"
#		command = "/usr/bin/time sleep 1 2>&1"

	puts command

#	time_output = `#{command}`.chomp
#	time_output = `#{command}`.split.last
	output = `#{command}`

#	       82.36 real        74.28 user         7.33 sys

#	irb(main):011:0> s
#	=> "       82.36 real        74.28 user         7.33 sys"
#	irb(main):012:0> s.split
#	=> ["82.36", "real", "74.28", "user", "7.33", "sys"]

#	may contain ...
#	Warning: Number of descriptions overridden to 10, number of alignments overridden to 10.

	
	puts output
	puts

	time_parts = output.split("\n").last.split("\s")
	
	csv = [iteration,fasta,evalue,num_threads,strand,num_descriptions,num_alignments,
		max_target_seqs,word_size,ungapped,outfmt,time_parts[0],time_parts[2],time_parts[4]]

	puts csv.join(',')
	puts

	c << csv

end	#	%w( 1e-1 1e-10 1e-100 ).each do |evalue|
end	#	[0,5].each do |outfmt|
end	#	[true,false].each do |ungapped|
end	#	).each do |word_size|
end	#	).each do |max_target_seqs|
end	#	).each do |num_alignments|
end	#	).each do |num_descriptions|
end	#	).each do |strand|
end	#	).each do |num_threads|
end	#	).each do |fasta|
end	#	( 1..10 ).each do |iteration|
c.close
__END__

jakewendt@fxdgroup-169-229-196-225 : blast_testing 619> time blastn -db /Users/jakewendt/blast_indexes/nt -num_threads 1 -query blast_input_3_reads_201.fasta -use_index false -evalue 1e-40 -strand both -num_descriptions 10 -num_alignments 10 -max_target_seqs 10 -ungapped

   Default = `-'
 -evalue <Real>
   Expectation value (E) threshold for saving hits 
   Default = `10'
 -word_size <Integer, >=4>
   Word size for wordfinder algorithm (length of best perfect match)


makembindex volsize 1 ... 660 ~68MB files
makembindex volsize 10 ... 1606 ... ~83MB files ?????




Simply having an index present triggers something which results in failures.  Yay!
Unless explicitly add "-use_index false" ?

Assertion failed: (0), function LocateIndex, file /Users/jakewendt/ncbi_cxx--7_0_0/src/algo/blast/api/blast_dbindex.cpp, line 195.

Perhaps because it is incomplete






/usr/bin/time blastn -db /Users/jakewendt/blast_indexes/nt -num_threads 1 -query blast_input_1_read_201a.fasta -evalue 1e-10 -use_index false -strand both -num_descriptions 10 -num_alignments 100 -max_target_seqs 10 -outfmt 0 -out /dev/null 2>&1
Warning: Number of descriptions overridden to 10, number of alignments overridden to 10.

num_descriptions, num_alignments and max_target_seqs appear to be related.
If one is out of sync, it will be adjusted.

See ...
	src/algo/blast/blastinput/blast_input_aux.cpp


