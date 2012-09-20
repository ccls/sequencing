#!/usr/bin/env ruby

# rins.rb is a tool to Rapidly Identify Nonhuman Sequences
#	(This is a CCLS modified version of the Stanford rins.pl by Jakeb)
# 
# script is written by Kun Qu and Aparna Bhaduri in Stanford Dermatology
#
#	rubified by jake
#

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'

class String
	def execute
		puts "Executing ..."
		puts self
		#	status is true or false
		status = system(self);
		#	$? is the pid and the return code
		puts ( "\n#{self} failed with #{$?}\n" ) unless ( status );
	end
end

#
#	This could be useful, but would have to also
#	specify the working directory.
#
start_step      = 0	

o = {
	:config_filename => 'config.yml',
	:output_filename => 'results.txt',
	:link_sample_fa_files => false,
	:chop_read_length => 25,
	:minIdentity => 80,
	:compress_ratio_thrd => 0.5,
	:iteration => 2,
	:bowtie_threads => 6,
	:bowtie_mismatch => 3,
	:paired_fragment_length => 300,
	:min_contig_length => 300,
	:trinity_threads => 6,
	:blastn_evalue_thrd => 0.05,
	:similarity_thrd => 0.8,
	:mailto => '',
	:die_on_failed_file_check => false,
	:files => {}
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	opts.banner = "Usage: #{$0} [options]\n" <<
		"Run the RINS pipeline to identify nonhuman sequences.\n" <<
		"Command Line Example: rins.pl -c config.txt -o output.txt\n" <<
		"	-h, --help    Displays this information\n" <<
		"	-c, --config  Config Filename (default: config.txt)\n" <<
		"	-o, --output  Output Filename (default: results.txt)\n" #<<
#		"	-s, --start   Start Step (default: 0)\n"

	# Define the options, and what they do

#	#	How to force this to be an integer?  Just cast it?
#	opts.on( '-s', '--start INTEGER', 'Start at step #' ) do |s|
#		start_step = s
#	end

	opts.on( '-c', '--config FILENAME', 'Processing options in ...' ) do |s|
		o[:config_filename] = s
	end

	opts.on( '-o', '--output FILENAME', 'Final output written to ...' ) do |s|
		o[:output_filename] = s
	end

end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options. What's left is the list of files to resize.
optparse.parse!

system "date";

config = YAML::load( ERB.new( IO.read( File.join( o[:config_filename] ) ) ).result)
o.update( config )

class RINS

	attr_accessor :options

	def initialize(options={})
		self.options = options
	end

	def method_missing(symb,*args,&block)
		if options.has_key? symb.to_s.to_sym
			options[symb.to_s.to_sym]
		else
			super
		end
	end

	def file_check( filename, empty_size = 0 )
	
	#	STDERR is not logged when using " | tee -a log"
	
		msg = '';
		unless( File.exists?(filename) )
			msg = "#{filename} not created";
			if( die_on_failed_file_check )
				raise msg;
			else
				puts msg
			end
		end
		#	Need to check existance too because if not dying on failure
		#	could have passed above check.  Then '-s $filename' is '' raising ...
		#	Use of uninitialized value in numeric le (<=) at ...
		if( ( File.exists?(filename) ) && ( File.size(filename) <= empty_size ) )
			msg = "#{filename} empty ( <= #{empty_size} )";
			if( die_on_failed_file_check )
				raise msg
			else
				puts msg
			end
		end
	end

	def step1
		raise "File format can either be fastq or fasta" unless( 
			['fasta','fastq'].include?(file_format) )
		
		outdir = "#{Time.now.strftime("%Y%m%d%H%M%S")}.outdir"
		FileUtils.mkdir outdir
		FileUtils.chdir outdir
		
		puts "step 1 change fastq files to fasta files"
		if( file_format == "fastq")
			files.each_pair do |k,v|
				"fastq2fasta.pl #{v} #{k}lane.fa".execute
				file_check( "#{k}lane.fa" );
			end
		else #if ($file_format eq "fasta") {
			if( link_sample_fa_files )
				puts "already fasta format, linking fasta files instead"
				files.each_pair do |k,v|
					FileUtils.ln_s(v,"#{k}lane.fa")
					file_check( "#{k}lane.fa" );
				end
			else
				puts "already fasta format, copy fasta files instead"
				files.each_pair do |k,v|
					FileUtils.cp(v,"#{k}lane.fa")
					file_check( "#{k}lane.fa" );
				end
			end
		end
	end

	def step2
		#
		#	TODO for some reason, blat doesn't work on the chopped????
		#
		#	blat34(2007) works, blat35(2012) does not???  Substantial diffs somewhere.
		#
		puts "step 2 chop reads"
		files.each_pair do |k,v|
			"chopreads.pl #{k}lane.fa chopped_#{k}lane.fa #{chop_read_length}".execute
			file_check( "chopped_#{k}lane.fa" )
		end
	end

	def step3
		puts "step 3 blat chopped reads"
		files.each_pair do |k,v|
			"blat #{blat_reference} -minIdentity=#{minIdentity} chopped_#{k}lane.fa chopped_#{k}lane.psl".execute
			file_check( "chopped_#{k}lane.psl", 427 )
		end
	end

	def step4
		puts "step 4 find blat out candidate reads"
		command = "blatoutcandidate.pl "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "chopped_#{k}lane.psl " }
		files.keys.sort.each{|k| command << "#{k}lane.fa " }
		command.execute
		files.each_pair { |k,v| file_check( "blat_out_candidate_#{k}lane.fa" ) }
	end

	def step5
		puts "step 5 compress raw reads"
		files.each_pair do |k,v|
			"compress.pl blat_out_candidate_#{k}lane.fa #{compress_ratio_thrd} > compress_#{k}lane.names".execute
			file_check( "compress_#{k}lane.names" )
		end
	end

	def step6
		puts "step 6 pull reads from blat_out_candidate fasta files"
		command = "pull_reads_fasta.pl "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "compress_#{k}lane.names " }
		files.keys.sort.each{|k| command << "blat_out_candidate_#{k}lane.fa " }
		files.keys.sort.each{|k| command << "compress_#{k}lane.fa " }
		command.execute
		files.each_pair { |k,v| file_check( "compress_#{k}lane.fa" ) }
	end

	def step7
		puts "step 7 align compressed reads to human genome reference using bowtie"
		files.each_pair do |k,v|
			"bowtie -n #{bowtie_mismatch} -p #{bowtie_threads} -f -S #{bowtie_index_human} compress_#{k}lane.fa compress_#{k}lane.sam".execute
			file_check( "compress_#{k}lane.sam" )
			"sam2names.pl compress_#{k}lane.sam bowtie_#{k}lane.names".execute
			file_check( "bowtie_#{k}lane.names" )
		end
		command = "pull_reads_fasta.pl "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "bowtie_#{k}lane.names " }
		files.keys.sort.each{|k| command << "compress_#{k}lane.fa " }
		files.keys.sort.each{|k| command << "bowtie_#{k}lane.fa " }
		command.execute
		files.each_pair { |k,v| file_check( "bowtie_#{k}lane.fa" ) }

		command = "candidate_non_human.pl "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "bowtie_#{k}lane.names " }
		command.execute
		file_check( "candidate_non_human.txt" )
	end

	def trinity_process
	
		#	--paired_fragment_length changed to --group_pairs_distance
		#	--run_butterfly no longer needed as is default
		#	--JM 1G now required
		
		#	changed --JM from 10G to just 1G and it was SO MUCH FASTER! NO DISK CHATTER!
		#		I don't have 10G of memory so it was probably using the disk as virtual memory.
		#		BIG MISTAKE.  Complete execution time is 23 minutes!  Still wrong though.
		
		#	--compatible_path_extension (for butterfly)
		# > java -Xmx20G -Xms1G -jar /Users/jakewendt/rins/trinity/Butterfly/Butterfly.jar -N 4817 -L 300 -F 300 -C /Users/jakewendt/rins_test/201209060909/trinity_output/chrysalis/RawComps.1/comp3552 --compatible_path_extension --stderr --max_number_of_paths_per_node=10 --path_reinforcement_distance=75 --triplet-lock
		#Exception in thread "main" java.lang.NullPointerException
		#	at gnu.getopt.Getopt.checkLongOption(Getopt.java:869)
		#	at gnu.getopt.Getopt.getopt(Getopt.java:1119)
		#	at TransAssembly_allProbPaths.main(TransAssembly_allProbPaths.java:193)
		#
		#	The compatible_path_extension option is not in TransAssembly_allProbPaths.java (it is the DEFAULT)
		#	It was last used in trinityrnaseq-r20110519
		#		trinityrnaseq-r20110519/Butterfly/src/src/TransAssembly_allProbPaths.java
		#
		#
		#	BEGIN DUPLICATE CODE
		#
		print "de novo assembly using Trinity\n";
		command = "Trinity.pl --seqType fa --group_pairs_distance #{paired_fragment_length} --min_contig_length #{min_contig_length} --output trinity_output --CPU #{trinity_threads} --bfly_opts \"--stderr\" --JM 1G "
		files.each_pair { |k,v| command << "--#{k} bowtie_#{k}lane.fa " }
		command.execute
		file_check( 'trinity_output/Trinity.fasta' )
		FileUtils.cp("trinity_output/Trinity.fasta","Trinity.fasta")
		#
		#	This script just joins the sequence on a single line
		#	rather that having it in 60 character line segments.
		#	I don't know why this is needed, but nevertheless ...
		#
		"modify_trinity_output.pl Trinity.fasta".execute
		puts "blastn trinity output against human genome"
		"blastn -query=Trinity.fasta -db=#{blastn_index_human} -evalue #{blastn_evalue_thrd} -outfmt 6 > human_contig.txt".execute
		file_check( 'human_contig.txt' )	#	NOTE  don't know how big an "empty" one is
		puts "clean up blastn outputs"
		"blastn_cleanup.pl human_contig.txt Trinity.fasta clean_blastn.fa #{similarity_thrd}".execute
		file_check( 'clean_blastn.fa' );	#	NOTE  don't know how big an "empty" one is
		FileUtils.rm_r("trinity_output")
		#
		#	END DUPLICATE CODE
		#
	end

	def step8
		nth_iteration = 1
		puts "step 8 #{nth_iteration} iteration"
		trinity_process
		for nth_iteration in 2..iteration
			puts "step 8 #{nth_iteration} iteration"
			puts "blat chopped reads"
			files.each_pair do |k,v|
				"blat clean_blastn.fa -minIdentity=95 #{k}lane.fa #{k}lane.iteration.psl".execute
				file_check( "#{k}lane.iteration.psl", 427 )
			end
			puts "find blat out candidate reads"
		
			command = "blatoutcandidate.pl "
			files.keys.sort.each{|k| command << "#{k}lane.iteration.psl " }
			files.keys.sort.each{|k| command << "#{k}lane.fa " }
			command.execute
			files.each_pair do |k,v|
				file_check( "blat_out_candidate_#{k}lane.fa" )
				FileUtils.cp("blat_out_candidate_#{k}lane.fa","iteration_#{k}lane.fa")
			end
			trinity_process
		end
	end

	def step9
		puts "step 9 detect species of non human sequences"
		puts "blastn trinity output against non-human genome"
		FileUtils.cp("clean_blastn.fa","non_human_contig.fa")
		"blastn -query=non_human_contig.fa -db=#{blastn_index_non_human} -evalue #{blastn_evalue_thrd} -outfmt 6 > non_human_contig_blastn.txt".execute
		file_check( 'non_human_contig_blastn.txt' );	#	NOTE  don't know how big an "empty" one is
		
		files.each_pair do |k,v|
			"blat non_human_contig.fa -minIdentity=98 iteration_#{k}lane.fa #{k}lane.psl".execute
			file_check( "#{k}lane.psl", 427 )
		end
		
		puts "write results"
		command = "write_result.pl non_human_contig.fa non_human_contig_blastn.txt "
		files.each_pair { |k,v| command << "#{k}lane.psl " }
		command << output_filename
		command.execute
	end

	def run
		step1
		step2
		step3
		step4
		step5
		step6
		step7
		step8
		step9
	end

end

app = RINS.new(o)
app.run



#	TODO	insert some sort of check

#&mailme();

system "date";

######################################################################

__END__


sub mailme {
	return if not defined $mailto;
	my $text = "RINS analysis completed";
	print "Attempting to mail $mailto the complete notice\n";
#	system "echo '$text' | mail -s 'RINS complete notice' $mailto";
}
