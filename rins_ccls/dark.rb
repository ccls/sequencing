#!/usr/bin/env ruby

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'
require 'ostruct'

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'
require 'file_format_detector'

puts "Running #{$0}"
o = {
	:config_filename          => 'config.yml',
	:output_filename          => 'results.txt',
	:output_suffix            => 'dark',
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
EOB

	# Define the options, and what they do

	opts.on( '--output_suffix STRING', 
		"Output directory suffix (default '#{o[:output_suffix]}')" ) do |s|
		o[:output_suffix] = s
	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end

	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
puts <<EOB
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

config = YAML::load( ERB.new( IO.read( File.join( o[:config_filename] ) ) ).result)
o.update( config )


#
#
#		Now we begin ....
#
#
system "date";

class Darkness < CclsSequencer

	def bowtie_non_human
		outbase = ''	#	outside to save last value for pull_reads_fasta
		files.each_pair do |k,v|	#	left,right
			outbase = 'lane'
			[bowtie_human_indexes].flatten.each do |bowtie_human_index|
				prevbase = String.new(outbase)	#	DO NOT JUST SET = AS WILL NOT CREATE NEW STRING!
				name = File.basename(bowtie_human_index)
				outbase << "_not_#{name}"
		
				command = "bowtie -n 3 -p 6 -f -S " <<
					"#{bowtie_human_index} #{k}#{prevbase}.fa " <<
					"#{k}#{outbase}.sam --un #{k}#{outbase}.fa"
				#
				#	add the --un for use with the next call
				#	I think that only the last .sam file will be needed to sync
				#
				command.execute
				"#{k}#{outbase}.sam".file_check(die_on_failed_file_check)
		
			end	#	%w( hg18 hg19 ).each do |hg|
			
			#
			#	sam2names only selects those names with a "*" in the third column.
			#	This is the column that contains the match name.  
			#	* means didn't match.
			#	Therefore, the output file is the sequences that did not match
			#	Not exactly the best filename-function.
			#
			"sam2names.rb #{k}#{outbase}.sam #{k}#{outbase}.names".execute
			"#{k}#{outbase}.names".file_check(die_on_failed_file_check)
		
		end	#	o.files.each_pair do |k,v|	#	left,right
		
		pull_reads_from_fastas(
			files.keys.sort.collect{|k| "#{k}#{outbase}.names" },
			files.keys.sort.collect{|k| "#{k}lane.fa" },
			files.keys.sort.collect{|k| "#{k}lane_bowtie_non_human.fa" })

		puts "de novo assembly using Trinity"
		command = "Trinity.pl --seqType fa " <<
			"--group_pairs_distance #{paired_fragment_length} " <<
			"--min_contig_length #{min_contig_length} " <<
			"--output trinity_output " <<
			"--CPU #{trinity_threads} " <<
			"--bfly_opts \"--stderr\" --JM 1G "

#			"--output trinity_output_#{nth_iteration} " <<
#
#	This is the only place where we NEED left and right or single.
#	Trinity expects the options.
#
		files.each_pair { |k,v| command << "--#{k} #{k}lane_bowtie_non_human.fa " }
		command.execute
		"trinity_output/Trinity.fasta".file_check(die_on_failed_file_check)
		FileUtils.cp("trinity_output/Trinity.fasta","trinity_non_human.fasta")
	end

end


#########

darkness = Darkness.new(o)
darkness.prepare_output_dir_and_log_file
darkness.file_format_check_and_conversion
darkness.bowtie_non_human
darkness.blastn_non_human("trinity_non_human.fasta")
darkness.wrap_things_up








__END__


	def align_compressed_reads_to_human_genome_reference_using_bowtie
		puts "step 7 align compressed reads to human genome reference using bowtie"
		files.each_pair do |k,v|
			#	bowtie's verbose is RIDICULOUS!
			#	It prints WAY too much and adds WAY too much time.
			#				"--verbose "<<
			command = "bowtie -n #{bowtie_mismatch} -p #{bowtie_threads} -f " <<
				"-S #{bowtie_index_human} compress_#{k}lane.fa compress_#{k}lane.sam"
			command.execute
			file_check( "compress_#{k}lane.sam" )	#	the reads that DIDN'T align?	NO

			"sam2names.rb compress_#{k}lane.sam bowtie_#{k}lane.names".execute
			file_check( "bowtie_#{k}lane.names" )
		end

		command = "pull_reads_fasta.rb "
		#	files is a hash and the keys are not guaranteed to be sorted
		#	sort alphabetically and left is first, right is last (conveniently)
		files.keys.sort.each{|k| command << "bowtie_#{k}lane.names " }  #	input
		files.keys.sort.each{|k| command << "compress_#{k}lane.fa " }   #	input
		files.keys.sort.each{|k| command << "bowtie_#{k}lane.fa " }     #	output
		command.execute
		files.each_pair { |k,v| file_check( "bowtie_#{k}lane.fa" ) }    #	non human?

#
#	This script has fixed input of chopped_leftlane.psl (and right or single)
#	BAD. BAD. BAD.	TODO
#	This is only informative and nothing uses the output
#	so could be commented out.
#
#
#	TODO Replaced with ruby version, but still in development
#
#
#		command = "candidate_non_human.rb "
#		#	files is a hash and the keys are not guaranteed to be sorted
#		#	sort alphabetically and left is first, right is last (conveniently)
#		files.keys.sort.each{|k| command << "bowtie_#{k}lane.names " }
#		command.execute
#		file_check( "candidate_non_human.txt" )
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
#
#	This is the only place where we NEED left and right or single.
#	Trinity expects the options.
#
		files.each_pair { |k,v| command << "--#{k} bowtie_#{k}lane.fa " }
		command.execute
		file_check(  "trinity_output_#{nth_iteration}/Trinity.fasta" )
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
		
			command = "blatoutcandidate.rb "
			files.keys.sort.each{|k| command << "#{k}lane.iteration.psl " }
			files.keys.sort.each{|k| command << "#{k}lane.fa " }	#	raw reads input
			command.execute
#
#	FYI This overwrites first blat_out_candidate files.
#
			files.each_pair do |k,v|
				#	
				#	raw reads with names in the psl files.
				#	
				file_check( "blat_out_candidate_#{k}lane.fa" )
#	why copy and not just move?
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
#	TODO now that I'm using my write_result.rb, should I 
#		add the add_descriptions_to_results functionality?
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

		puts "parsing write results' results and selecting just first match"
		command = "just_first_contig_match_results.rb " <<
			"-i #{output_filename}.with_descriptions " <<
			"-o #{output_filename}.with_descriptions.just_first_match"
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

	def wrap_things_up
		puts "Finished at ..."
		system("date")
		STDOUT.reopen(original_stdout)
		STDERR.reopen(original_stderr)
		puts "All done."
		system("date")
	end

	def run
		prepare_output_dir_and_log_file
		file_format_check_and_conversion
#		if mode == 1	
			chop_reads
			blat_chopped_reads
			blat_out_candidate_reads
			compress_raw_reads
			pull_reads_from_blat_out_candidates
#		elsif mode == 2
#			files.each_pair { |k,v| FileUtils.ln_s("#{k}lane.fa","compress_#{k}lane.fa") }
#		end
		align_compressed_reads_to_human_genome_reference_using_bowtie
		step8
		detect_species_of_non_human_sequences
		detect_unknown_sequences
		wrap_things_up
	end

end

app = RINS.new(o)
app.run

