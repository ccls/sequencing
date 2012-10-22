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

#
#
#	TODO
#
#		Better file names and make scripts follow some conventions.
#		If I'm gonna continue with this, I may have to re-write all
#		of this.
#
#

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'
require 'file_format_detector'

#
#	This could be useful, but would have to also
#	specify the working directory.
#
#	start_step      = 0	
#

#
#	Default values
#
#	:output_suffix            => '',
o = {
	:config_filename          => 'config.yml',
	:output_filename          => 'results.txt',
#	:mode                     => 1,
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

Run the RINS pipeline to identify nonhuman sequences.
Command Line Example: rins.pl -c config.txt -o output.txt
EOB

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

#	opts.on( '-m', '--mode INTEGER', 
#		"Development mode (default #{o[:mode]})" ) do |s|
#		o[:mode] = s.to_i
#	end

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

#	--mode INTEGER
#
#		By default this is 1, which is 'normal' RINS processing.
#		If set to 2, this will skip the chopping and initial blat calls.
#
	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
puts <<EOB

------------------

command line

	--output_suffix STRING

		By default this is blank and the working directory is just a date/time stamp.
		If given, this string will by joined to the date/time stamp with 
		a '.' between.

		Example:

			default
				working dir -> 20121015091927.outdir

			--output_suffix skip_blat
				working dir -> 20121015091927.outdir.skip_blat


------------------

config file (yml format does not allow TABS) ...

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
	
	#	MinIdentity match for the initial blat candidate call
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

	#	Used with the blastn_cleanup script to determine match quality.
	#	Values seem to range from 0 to 2 with 0 being very unlikely
	#	and 2 being very likely.
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

class RINS < CclsSequencer

	#	just an alias really
	def blat_references
		blat_reference
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
				"chopped_#{k}lane.fa".file_check(die_on_failed_file_check)
			else
				"chopreads.pl #{k}lane.fa chopped_#{k}lane.fa #{chop_read_length}".execute
				"chopped_#{k}lane.fa".file_check(die_on_failed_file_check)
			end
		end
	end

	def blat_chopped_reads
		puts "step 3 blat chopped reads to KNOWN NON-HUMAN REFERENCES"
		blat_refs = [blat_references].flatten
		files.each_pair do |k,v|
			blat_refs.each do |blat_ref|
				basename = File.basename(blat_ref)
				command = "blat #{blat_ref} -dots=1000 " <<
					"-minIdentity=#{minIdentity} " <<
					"chopped_#{k}lane.fa " <<
					"chopped_#{k}lane_#{basename}.psl"
				command.execute
#
#	"chopped_#{k}lane_#{basename}.psl" contains the candidates from the 
#		blat_reference.  In these cases, they should be non-human matches.
#
				"chopped_#{k}lane_#{basename}.psl".file_check(die_on_failed_file_check,427)
			end

			if blat_refs.length > 1
				puts "Merging chopped_#{k}lane psl files."
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
			else
				FileUtils.ln_s("chopped_#{k}lane_#{File.basename(blat_refs[0])}.psl", 
					"chopped_#{k}lane.psl")
			end
		end
	end

#
#	blatoutcandidate.pl loops over the lines in the given psl files and
#		gets the 10ths column ... @HWI-ST281_0133:3:1:6254:2049#0/1
#		removes the trailing "/0" or "/1" leaving @HWI-ST281_0133:3:1:6254:2049#0
#		and adds it to the hash
#
#	(our sequence names don't match this /0 or /1 so I wrote our own ruby version)
#
#	Then it loops of the first fasta file 
#		(left hopefully as it writes to blat_out_candidate_leftlane.fa)
#		If the sequence line name is in the hash,
#			prints it and the sequence to the output fasta file
#	Then same thing for the second fasta file
#		(right hopefully as it writes to blat_out_candidate_rightlane.fa)
#
#	Interesting that both psl files are actually used for each fasta file.
#

	def find_blat_out_candidate_reads
		puts "step 4 find blat out candidate reads"
		blat_out_candidate_reads(
			files.keys.sort.collect{|k| "chopped_#{k}lane.psl " },
			files.keys.sort.collect{|k| "#{k}lane.fa " },
			files.keys.sort.collect{|k| "04_blat_out_candidate_#{k}lane.fa" })

#		command = "blatoutcandidate.rb "
#		#	files is a hash and the keys are not guaranteed to be sorted
#		#	sort alphabetically and left is first, right is last (conveniently)
#		files.keys.sort.each{|k| command << "chopped_#{k}lane.psl " } #	NON-HUMAN matches
#		files.keys.sort.each{|k| command << "#{k}lane.fa " } #	raw reads input
#		command.execute
##
##	blatoutcandidate.pl ALWAYS creates ... blat_out_candidate_#{k}lane.fa
##	I REALLY don't like that.  So much inconsistancy. Will overwrite existing.
##
##	TODO wrote my own version of blatoutcandidate so could change this
##
#		files.each_pair { |k,v| 
#			#	
#			#	raw reads with names in the psl files.
#			#	
#			"blat_out_candidate_#{k}lane.fa".file_check(die_on_failed_file_check)
#			FileUtils.mv( "blat_out_candidate_#{k}lane.fa",
#				"04_blat_out_candidate_#{k}lane.fa" )	#	NON-HUMAN matches 
#		}
	end

	def compress_raw_reads
		puts "step 5 compress raw reads"
		files.each_pair do |k,v|
			command = "compress.pl 04_blat_out_candidate_#{k}lane.fa #{compress_ratio_thrd} " <<
				"> compress_#{k}lane.names"
			command.execute
			"compress_#{k}lane.names".file_check(die_on_failed_file_check) # NON-HUMAN matches
		end
	end

	def pull_reads_from_blat_out_candidates
		puts "step 6 pull reads from blat_out_candidate fasta files"
		pull_reads_from_fastas(
			files.keys.sort.collect{|k| "compress_#{k}lane.names" },
			files.keys.sort.collect{|k| "04_blat_out_candidate_#{k}lane.fa" },
			files.keys.sort.collect{|k| "compress_#{k}lane.fa" }) # NON-HUMAN matches
	end


#
#	The bowtie created SAM file includes all sequences with info regarding
#	its match or non-match.  The sam2names script selects only those
#	entries which contain a "*" in what I suppose is the appropriate column
#	that flags it as having not matched.  These unaligned sequence names
#	are then pulled by pull_reads_fasta.rb to create the new filtered
#	and paired fasta files.
#

	def align_compressed_reads_to_human_genome_reference_using_bowtie
		puts "step 7 align compressed reads to human genome reference using bowtie"
		files.each_pair do |k,v|
			#	bowtie's verbose is RIDICULOUS!
			#	It prints WAY too much and adds WAY too much time.
			#				"--verbose "<<
			command = "bowtie -n #{bowtie_mismatch} -p #{bowtie_threads} -f " <<
				"-S #{bowtie_index_human} compress_#{k}lane.fa compress_#{k}lane.sam"
			command.execute
			"compress_#{k}lane.sam".file_check(die_on_failed_file_check) #	the reads that DIDN'T align?	NO

			"sam2names.rb compress_#{k}lane.sam bowtie_#{k}lane.names".execute
			"bowtie_#{k}lane.names".file_check(die_on_failed_file_check)
		end

		pull_reads_from_fastas(
			files.keys.sort.collect{|k| "bowtie_#{k}lane.names" },
			files.keys.sort.collect{|k| "compress_#{k}lane.fa" },
			files.keys.sort.collect{|k| "bowtie_#{k}lane.fa" })

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
		puts "de novo assembly using Trinity"
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
		"trinity_output_#{nth_iteration}/Trinity.fasta".file_check(die_on_failed_file_check)
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
		'human_contig.txt'.file_check(die_on_failed_file_check)	#	NOTE  don't know how big an "empty" one is

		puts "clean up blastn outputs"
		command = "blastn_cleanup.rb human_contig.txt Trinity.fasta " <<
			"clean_blastn.fa #{similarity_thrd}"
		command.execute
		'clean_blastn.fa'.file_check(die_on_failed_file_check) #	NOTE  don't know how big an "empty" one is
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
				"#{k}lane.iteration.psl".file_check(die_on_failed_file_check,427)
			end
#			puts "find blat out candidate reads"
		
			blat_out_candidate_reads(
				files.keys.sort.collect{|k| "#{k}lane.iteration.psl " },
				files.keys.sort.collect{|k| "#{k}lane.fa " },
				files.keys.sort.collect{|k| "iteration_#{k}lane.fa" })

#			command = "blatoutcandidate.rb "
#			files.keys.sort.each{|k| command << "#{k}lane.iteration.psl " }
#			files.keys.sort.each{|k| command << "#{k}lane.fa " }	#	raw reads input
#			command.execute
##
##	FYI This overwrites first blat_out_candidate files.
##
#			files.each_pair do |k,v|
#				#	
#				#	raw reads with names in the psl files.
#				#	
#				"blat_out_candidate_#{k}lane.fa".file_check(die_on_failed_file_check)
##	why copy and not just move?
#				FileUtils.cp("blat_out_candidate_#{k}lane.fa","iteration_#{k}lane.fa")
#			end
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
		'non_human_contig_blastn.txt'.file_check(die_on_failed_file_check) #	NOTE  don't know how big an "empty" one is
		
		files.each_pair do |k,v|
			command = "blat non_human_contig.fa -minIdentity=98 " <<
				"iteration_#{k}lane.fa #{k}lane.psl"
			command.execute
			"#{k}lane.psl".file_check(die_on_failed_file_check,427)
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

end

app = RINS.new(o)
app.prepare_output_dir_and_log_file
app.file_format_check_and_conversion
#if mode == 2
#		app.files.each_pair { |k,v| FileUtils.ln_s("#{k}lane.fa","compress_#{k}lane.fa") }
#else #if mode == 1	
		app.chop_reads
		app.blat_chopped_reads
		app.find_blat_out_candidate_reads
		app.compress_raw_reads
		app.pull_reads_from_blat_out_candidates
#end
app.align_compressed_reads_to_human_genome_reference_using_bowtie
app.step8
app.detect_species_of_non_human_sequences
app.detect_unknown_sequences
app.wrap_things_up

