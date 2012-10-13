#!/usr/bin/env ruby

# rins.rb is a tool to Rapidly Identify Nonhuman Sequences
#	(This is a CCLS modified version of the Stanford rins.pl by Jakeb)
# 
# script is written by Kun Qu and Aparna Bhaduri in Stanford Dermatology
#
#	rubified by jake
#
#	NOTE THIS REQUIRES RUBY 1.9
#		for at least 1 reason... an array of symbols doesn't sort in 1.8.7
#			>> [:left,:right].sort
#			NoMethodError: undefined method `<=>' for :left:Symbol
#				from (irb):1:in `sort'
#				from (irb):1
#

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'file_format_detector'

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

#
#	Default values
#
#	:output_suffix            => '',
o = {
	:config_filename          => 'config.yml',
	:output_filename          => 'results.txt',
	:mode                     => 1,
	:link_sample_fa_files     => false,
	:pre_chopped              => false,
	:chop_read_length         => 25,
	:minIdentity              => 80,
	:compress_ratio_thrd      => 0.5,
	:iteration                => 2,
	:bowtie_threads           => 6,
	:bowtie_mismatch          => 3,
	:paired_fragment_length   => 300,
	:min_contig_length        => 300,
	:trinity_threads          => 6,
	:blastn_evalue_thrd       => 0.05,
	:similarity_thrd          => 0.8,
	:mailto                   => '',
	:die_on_failed_file_check => false,
	:files                    => {}
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB
Usage: #{File.basename($0)} [options]
   or: #{File.basename($0)} [options] | tee -a log &

Run the RINS pipeline to identify nonhuman sequences.
Command Line Example: rins.pl -c config.txt -o output.txt
EOB

#	opts.banner = "\nUsage: #{File.basename($0)} [options]\n\n" <<
#		"Run the RINS pipeline to identify nonhuman sequences.\n" <<
#		"Command Line Example: rins.pl -c config.txt -o output.txt\n\n" <<
#		"In the config file ...\n\n" <<
#		"	(:files is a hash)\n" <<
#		"	:files:\n" <<
#		"	  :left:  somefile1.fa\n" <<
#		"	  :right: somefile2.fa\n" <<
#		"		OR\n" <<
#		"	:files:\n" <<
#		"	  :single:  somefile.fa\n\n" <<
#		"	(:blat_reference CAN be an array)\n" <<
#		"	:blat_reference:\n" <<
#		"	  - somefile1.fa\n" <<
#		"	  - somefile2.fa\n" <<
#		"		OR\n" <<
#		"	:blat_reference:\n" <<
#		"	  - somefile.fa\n" <<
#		"		OR\n" <<
#		"	:blat_reference: somefile.fa\n\n" <<
#		"------------------\n\n"


	# Define the options, and what they do

#	#	How to force this to be an integer?  Just cast it?
#	opts.on( '-s', '--start INTEGER', 'Start at step #' ) do |s|
#		start_step = s
#	end

	opts.on( '-c', '--config FILENAME', 
		"Processing options in (default #{o[:config_filename]})" ) do |s|
		o[:config_filename] = s
	end

	opts.on( '-o', '--output FILENAME', 
		"Final output written to (default #{o[:output_filename]})" ) do |s|
		o[:output_filename] = s
	end

	opts.on( '-m', '--mode INTEGER', 
		"Development mode (default #{o[:mode]})" ) do |s|
		o[:mode] = s.to_i
	end

	opts.on( '--output_suffix STRING', 
		"Output directory suffix (default '#{o[:output_suffix]}')" ) do |s|
		o[:output_suffix] = s
	end

	# This displays the help screen, all programs are assumed to have this option.
	#	Add extra "\n" to last option for aesthetics.
#	opts.on( '-h', '--help', 'Display this screen' ) do
#		puts opts
#		exit
#	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end

	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
puts <<EOB

------------------

In the config file ...

	(:files IS a hash)
	:files:
	  :left:  somefile1.fa
	  :right: somefile2.fa
		OR
	:files:
	  :single:  somefile.fa

	(:blat_reference CAN be an array)
	:blat_reference:
	  - somefile1.fa
	  - somefile2.fa
		OR
	:blat_reference:
	  - somefile.fa
		OR
	:blat_reference: somefile.fa


	#	if the raw files are fasta, link them instead of copying them in
	:link_sample_fa_files: true


	#
	#	Link if pre_chopped, the chopped filename must be the same as the raw
	#	file name with "_chopped_(chop_read_length)" at end before extension.
	#
	#	testset_forkun_1.fa -> testset_forkun_1_chopped_25.fa
	#
	:pre_chopped: true


	# reference files
	#
	#	Blat requires the file suffix.
	#	The other references do not.
	#	There is a 4GB file size limit in blat.
	#	all_bacterial.fna is 4GB.
	#	virus_all.fasta is just a few hundred meg.
	#	Combining the 2 will raises errors.
	#	So, the blat_reference can be an array.
	#	blat will be run on each and the results combined.
	#
	:blat_reference:         
	  - /Volumes/cube/working/indexes/all_bacterial.fna
	  - /Volumes/cube/working/indexes/virus_all.fasta

	:bowtie_index_human:     /Users/jakewendt/rins/indexes/hg18
	:blastn_index_human:     /Users/jakewendt/rins/indexes/hg18
	:blastn_index_non_human: /Users/jakewendt/rins/indexes/virus


	#
	#	Additional options with their defaults include ...
	#
	#	The results output filename 
	#	(can be specified on the command line or in the config file)
	#		:output_filename: results.txt
	#
	#	Size of the chopped reads for the initial blat call
	:chop_read_length: 25
	
	#	MinIdentity match for the initial blat call
	:minIdentity: 80
	
	#
	:compress_ratio_thrd: 0.5
	
	#	The number of trinity iterations that are run
	:iteration: 2
	
	#
	:bowtie_threads: 6
	
	#
	:bowtie_mismatch: 3
	
	#
	:paired_fragment_length: 300
	
	#
	:min_contig_length: 300
	
	#
	:trinity_threads: 6
	
	#
	:blastn_evalue_thrd: 0.05

	#
	:similarity_thrd: 0.8

	#	Exit if an expected file is empty or missing
	:die_on_failed_file_check: false
	
EOB
	exit
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

	#	just an alias really
	def blat_references
		blat_reference
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

	def file_format_check_and_conversion


#
#	use file_format_detector
#
#		file_format = FileFormatDetector.new(files.first.value).format
#


		raise "File format can either be fastq or fasta" unless( 
			['fasta','fastq'].include?(file_format) )
		
		outdir = ["#{Time.now.strftime("%Y%m%d%H%M%S")}.outdir",options[:output_suffix]].compact.join('.')
		FileUtils.mkdir outdir
		FileUtils.chdir outdir
		
		puts "step 1 change fastq files to fasta files"
		files.each_pair do |k,v|
			if( file_format == "fastq")
				"fastq2fasta.pl #{v} #{k}lane.fa".execute
				file_check( "#{k}lane.fa" );
			else #if ($file_format eq "fasta") {
				if( link_sample_fa_files )
					puts "already fasta format, linking #{v} #{k}lane.fa instead"
					FileUtils.ln_s(v,"#{k}lane.fa")
					file_check( "#{k}lane.fa" );
				else
					puts "already fasta format, copying #{v} #{k}lane.fa instead"
					FileUtils.cp(v,"#{k}lane.fa")
					file_check( "#{k}lane.fa" );
				end
			end
		end
	end

	def chop_reads
		#
		#	TODO for some reason, blat doesn't work on the chopped????
		#
		#	blat34(2007) works, blat35(2012) does not???  Substantial diffs somewhere.
		#
		puts "step 2 chop reads"
		files.each_pair do |k,v|
			if( pre_chopped )
#
#	The raw files cannot be the same as the chopped files or
#	trinity crashes.  Trying naming convention.
#
				#	.fa or .fasta should both work
				chopped = v.gsub(/\.fa/,"_chopped_#{chop_read_length}.fa")
				puts "files are pre-chopped so linking #{chopped} chopped_#{k}lane.fa"

				FileUtils.ln_s(chopped,"chopped_#{k}lane.fa")
##				FileUtils.ln_s("#{k}lane.fa","chopped_#{k}lane.fa")
#puts "trinity having problems?  gonna copying chopped file"
#FileUtils.cp("#{k}lane.fa","chopped_#{k}lane.fa")
				file_check( "chopped_#{k}lane.fa" )
			else
				"chopreads.pl #{k}lane.fa chopped_#{k}lane.fa #{chop_read_length}".execute
				file_check( "chopped_#{k}lane.fa" )
			end
		end
	end

	def blat_chopped_reads
		puts "step 3 blat chopped reads"
		blat_refs = [blat_references].flatten
		files.each_pair do |k,v|
			blat_refs.each do |blat_ref|
				basename = File.basename(blat_ref)
				command = "blat #{blat_ref} -dots=1000 " <<
					"-minIdentity=#{minIdentity} " <<
					"chopped_#{k}lane.fa " <<
					"chopped_#{k}lane_#{basename}.psl"
				command.execute
				file_check( "chopped_#{k}lane_#{basename}.psl", 427 )
			end
	
			puts "Copying chopped_#{k}lane_#{File.basename(blat_refs[0])}.psl " <<
				"to chopped_#{k}lane.psl"
			FileUtils.cp("chopped_#{k}lane_#{File.basename(blat_refs[0])}.psl", 
				"chopped_#{k}lane.psl")
			( blat_refs - [blat_refs[0]] ).each do |blat_ref|
#
#	The blat header is 5 lines, not just 1
#
				command = "tail +6 chopped_#{k}lane_#{File.basename(blat_ref)}.psl " <<
					">> chopped_#{k}lane.psl"
				command.execute
			end
		end

	end

	def blat_out_candidate_reads
		puts "step 4 find blat out candidate reads"
		command = "blatoutcandidate.pl "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "chopped_#{k}lane.psl " }
		files.keys.sort.each{|k| command << "#{k}lane.fa " }
		command.execute
		files.each_pair { |k,v| file_check( "blat_out_candidate_#{k}lane.fa" ) }
	end

	def compress_raw_reads
		puts "step 5 compress raw reads"
		files.each_pair do |k,v|
			command = "compress.pl blat_out_candidate_#{k}lane.fa #{compress_ratio_thrd} " <<
				"> compress_#{k}lane.names"
			command.execute
			file_check( "compress_#{k}lane.names" )
		end
	end

	def pull_reads_from_blat_out_candidates
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

	def align_compressed_reads_to_human_genome_reference_using_bowtie
		puts "step 7 align compressed reads to human genome reference using bowtie"
		files.each_pair do |k,v|
#	bowtie's verbose is RIDICULOUS!
#	It prints way too much and adds way too much time.
#				"--verbose "<<
			command = "bowtie -n #{bowtie_mismatch} -p #{bowtie_threads} -f " <<
				"-S #{bowtie_index_human} compress_#{k}lane.fa compress_#{k}lane.sam"
			command.execute
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

	def trinity_process(nth_iteration)
	
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
		print "de novo assembly using Trinity\n";
		command = "Trinity.pl --seqType fa " <<
			"--group_pairs_distance #{paired_fragment_length} " <<
			"--min_contig_length #{min_contig_length} " <<
			"--output trinity_output_#{nth_iteration} " <<
			"--CPU #{trinity_threads} " <<
			"--bfly_opts \"--stderr\" --JM 1G "
		files.each_pair { |k,v| command << "--#{k} bowtie_#{k}lane.fa " }
		command.execute
		file_check( "trinity_output_#{nth_iteration}/Trinity.fasta" )
		FileUtils.cp("trinity_output_#{nth_iteration}/Trinity.fasta","Trinity.fasta")
		#
		#	This script just joins the sequence on a single line
		#	rather that having it in 60 character line segments.
		#	I don't know why this is needed, but nevertheless ...
		#
		#	modify_trinity_output.pl just puts the sequence on a single line
		#	This is kinda needed by blastn_cleanup.rb (not using the old perl version anymore)
		"modify_trinity_output.pl Trinity.fasta".execute
		puts "blastn trinity output against human genome"
		command = "blastn -query=Trinity.fasta -db=#{blastn_index_human} " <<
			"-evalue #{blastn_evalue_thrd} -outfmt 6 > human_contig.txt"
		command.execute
		file_check( 'human_contig.txt' )	#	NOTE  don't know how big an "empty" one is
		puts "clean up blastn outputs"
		command = "blastn_cleanup.rb human_contig.txt Trinity.fasta " <<
			"clean_blastn.fa #{similarity_thrd}"
		command.execute
		file_check( 'clean_blastn.fa' );	#	NOTE  don't know how big an "empty" one is
#		FileUtils.rm_r("trinity_output")
#
#	Why move it?  Why not just process it there?
#
#		FileUtils.mv("trinity_output","trinity_output_#{nth_iteration}")
	end

	def step8
		nth_iteration = 1
		puts "step 8 #{nth_iteration} iteration"
		trinity_process(nth_iteration)
		for nth_iteration in 2..iteration
			puts "step 8 #{nth_iteration} iteration"
			puts "blat chopped reads"
			files.each_pair do |k,v|
				command = "blat clean_blastn.fa -minIdentity=95 #{k}lane.fa " <<
					"#{k}lane.iteration.psl"
				command.execute
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
			trinity_process(nth_iteration)
		end
	end

	def detect_species_of_non_human_sequences
		puts "step 9 detect species of non human sequences"
		puts "blastn trinity output against non-human genome"
		FileUtils.cp("clean_blastn.fa","non_human_contig.fa")
		command = "blastn -query=non_human_contig.fa -db=#{blastn_index_non_human} " <<
			"-evalue #{blastn_evalue_thrd} -outfmt 6 > non_human_contig_blastn.txt"
		command.execute
		file_check( 'non_human_contig_blastn.txt' );	#	NOTE  don't know how big an "empty" one is
		
		files.each_pair do |k,v|
			command = "blat non_human_contig.fa -minIdentity=98 " <<
				"iteration_#{k}lane.fa #{k}lane.psl"
			command.execute
			file_check( "#{k}lane.psl", 427 )
		end
		
		puts "write results"
#		command = "write_result.pl non_human_contig.fa non_human_contig_blastn.txt "
#
#	Using my ruby write result now
#
		command = "write_result.rb non_human_contig.fa " <<
			"non_human_contig_blastn.txt "
		files.each_pair { |k,v| command << "#{k}lane.psl " }
		command << output_filename
		command.execute

		puts "parsing write results' results and adding a description"
#		command = "add_descriptions_to_results.rb -i #{output_filename} -o #{output_filename}.with_descriptions -d /Volumes/cube/working/indexes/all_bacterial_and_viral"
#		command = "add_descriptions_to_results.rb -i #{output_filename} -o #{output_filename}.with_descriptions"
#
#	the database needs to be the same as the one blastn used
#	otherwise, the sequence may not be in it.
#
		command = "add_descriptions_to_results.rb " <<
			"-i #{output_filename} " <<
			"-o #{output_filename}.with_descriptions " <<
			"-d #{blastn_index_non_human}"
		command.execute
	end

	def detect_unknown_sequences
		puts "step X detect the unknowns (this is experimental)"
#	make sure sequences are on single line
		command = "blastn_cleanup.rb " <<
			"non_human_contig_blastn.txt " <<
			"non_human_contig.fa " <<
			"unknown_sequences.fa 0.5"
		command.execute
	end

	def run
		file_format_check_and_conversion
		if mode == 1	
			chop_reads
			blat_chopped_reads
			blat_out_candidate_reads
			compress_raw_reads
			pull_reads_from_blat_out_candidates
		elsif mode == 2
			files.each_pair { |k,v| FileUtils.ln_s("#{k}lane.fa","compress_#{k}lane.fa") }
		end
		align_compressed_reads_to_human_genome_reference_using_bowtie
		step8
		detect_species_of_non_human_sequences
		detect_unknown_sequences
		puts "All done."
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
