#!/usr/bin/env ruby

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'
require 'ostruct'

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'

puts "Running #{$0}"
o = {
	:config_filename          => 'config.yml',
	:output_filename          => 'results.txt',
	:output_suffix            => 'dark',
	:bowtie_version           => 1,
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

class Darkness < CclsSequencer

#
#	This was designed to work with singles and pairs.
#	Now it really kinda requires pairs.
#

	def bowtie_non_human
		outbase = ''	#	outside to save last value for pull_reads_fasta

		outbase = 'raw'
		files.keys.sort.each_with_index{|k,i| 
			FileUtils.ln_s(files[k],"raw.#{i+1}.#{file_format}") }
		"raw.1.#{file_format}".file_check(die_on_failed_file_check)
		"raw.2.#{file_format}".file_check(die_on_failed_file_check)

		[bowtie_human_indexes].flatten.each_with_index do |bowtie_human_index,i|
			prevbase = String.new(outbase)	#	DO NOT JUST SET = AS WILL NOT CREATE NEW STRING!
			name = File.basename(bowtie_human_index)
			outbase << "_not" if i == 0
			outbase << "_#{name}"

			command = "bowtie2 -p 4 -N 1 " <<
				"#{(file_format == 'fastq')? '-q' : '-f'} " <<
				"-x #{bowtie_human_index} " <<
				"-1 #{prevbase}.1.#{file_format} " <<
				"-2 #{prevbase}.2.#{file_format} " <<
				"-S #{outbase}.sam " <<
				"--un-conc #{outbase}.#{file_format}"
			command.execute

			"#{outbase}.sam".file_check(die_on_failed_file_check)
			"#{outbase}.1.#{file_format}".file_check(die_on_failed_file_check)
			"#{outbase}.2.#{file_format}".file_check(die_on_failed_file_check)
	
		end	#	%w( hg18 hg19 ).each do |hg|
	
		
#		pull_reads_from_fastas(
#			files.keys.sort.collect{|k| "#{k}#{outbase}.names" },
#			files.keys.sort.collect{|k| "#{k}lane.fa" },
#			files.keys.sort.collect{|k| "#{k}lane_bowtie_non_human.fasta" })
#
#	link last files as 
#	files.keys.sort.collect{|k| "#{k}lane_bowtie_non_human.fasta" })

#		files.keys.sort.each_with_index{|k,i| 
#			FileUtils.ln_s("#{outbase}.#{i+1}.#{file_format}","#{k}lane_bowtie_non_human.#{file_format}") }

		puts "de novo assembly using Trinity"
#		command = "Trinity.pl --seqType fa " <<
		command = "Trinity.pl --seqType #{(file_format == 'fastq')? 'fq' : 'fa'} " <<
			"--group_pairs_distance #{paired_fragment_length} " <<
			"--min_contig_length #{min_contig_length} " <<
			"--output trinity_output " <<
			"--CPU #{trinity_threads} " <<
			"--bfly_opts \"--stderr\" --JM 1G "
		files.keys.sort.each_with_index{|k,i| 
			command << "--#{k} #{outbase}.#{i+1}.#{file_format} " }

#		files.each_pair { |k,v| command << "--#{k} #{k}lane_bowtie_non_human.#{file_format} " }


#
#	will fastq in produce fastq out with trinity????
#

		command.execute
		"trinity_output/Trinity.fasta".file_check(die_on_failed_file_check)
		FileUtils.cp("trinity_output/Trinity.fasta","trinity_non_human.fasta")

	end

end


#########

system "date";
darkness = Darkness.new(o)
darkness.prepare_output_dir_and_log_file


#darkness.file_format_check_and_conversion


darkness.bowtie_non_human
darkness.blastn_non_human("trinity_non_human.fasta")
darkness.wrap_things_up





__END__
