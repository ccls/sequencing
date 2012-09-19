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

def file_check( filename, empty_size = 0 )

#	STDERR is not logged when using " | tee -a log"

	msg = '';
	unless( File.exists?(filename) )
		msg = "#{filename} not created";
		if( $die_on_failed_file_check )
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

def do_this( command )
	puts "Executing ..."
	puts "#{command}"

#	backticks cache the output so you see nothing until its done.
#	@result would be that output that COULD be printed.
#	my @result = `$command`;	

#	$exit_status here is the same as the built in $?
#	so long as another system command hasn't been run.
#	my $exit_status = system($command);

	status = system(command);

	puts ( "\n#{command} failed with #{$?}\n" ) unless ( status );
end

config_filename = 'config.yml'
output_filename = 'results.txt'



#
#	This could be useful, but would have to also
#	specify the working directory.
#
start_step      = 0	




optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	opts.banner = "Usage: #{$0} [options]\n" <<
		"Run the RINS pipeline to identify nonhuman sequences.\n" <<
		"Command Line Example: rins.pl -c config.txt -o output.txt\n" <<
		"	-h, --help    Displays this information\n" <<
		"	-c, --config  Config Filename (default: config.txt)\n" <<
		"	-o, --output  Output Filename (default: results.txt)\n" <<
		"	-s, --start   Start Step (default: 0)\n"

	# Define the options, and what they do

	#	How to force this to be an integer?  Just cast it?
	opts.on( '-s', '--start INTEGER', 'Start at step #' ) do |s|
		start_step = s
	end

	opts.on( '-c', '--config FILENAME', 'Processing options in ...' ) do |s|
		config_filename = s
	end

	opts.on( '-o', '--output FILENAME', 'Final output written to ...' ) do |s|
		output_filename = s
	end

end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options. What's left is the list of files to resize.
optparse.parse!

system "date";

config = YAML::load( ERB.new( IO.read( File.join( config_filename ) ) ).result)

# Config values of files and files' configurations

file_format = config[:file_format]
#pair_end    = config[:pair_end]	#	irrelevant
link_sample_fa_files = config[:link_sample_fa_files] || false
#leftlane_filename	  = config[:leftlane_filename]	#	irrelevant
#rightlane_filename  = config[:rightlane_filename]	#	irrelevant
#singlelane_filename = config[:singlelane_filename]	#	irrelevant
chop_read_length = config[:chop_read_length] || 25
minIdentity = config[:minIdentity] || 80

blat_reference = config[:blat_reference]

compress_ratio_thrd = config[:compress_ratio_thrd] || 0.5
iteration = config[:iteration] || 2

# executables' configurations
blat_bin = config[:blat_bin] || 'blat'
bowtie_bin = config[:bowtie_bin] || 'bowtie'
bowtie_build_bin = config[:bowtie_build_bin] || 'bowtie-build'
bowtie_index_human = config[:bowtie_index_human]
bowtie_threads = config[:bowtie_threads] || 6
bowtie_mismatch = config[:bowtie_mismatch] || 3

trinity_script = config[:trinity_script] || 'Trinity.pl'
paired_fragment_length = config[:paired_fragment_length] || 300
min_contig_length = config[:min_contig_length] || 300
trinity_threads = config[:trinity_threads] || 6

blastn_bin = config[:blastn_bin] || 'blastn'
blastn_index_human = config[:blastn_index_human]
blastn_index_non_human = config[:blastn_index_non_human]
blastn_evalue_thrd = config[:blastn_evalue_thrd] || 0.05
similarity_thrd = config[:similarity_thrd] || 0.8


fastq2fasta_script = "fastq2fasta.pl";
chopreads_script = "chopreads.pl";
blat_out_candidate_script = "blatoutcandidate.pl";
compress_script = "compress.pl";
pull_reads_fasta_script = "pull_reads_fasta.pl";
sam2names_script = "sam2names.pl";
modify_trinity_output_script = "modify_trinity_output.pl";
blastn_cleanup_script = "blastn_cleanup.pl";
candidate_non_human_script = "candidate_non_human.pl";
write_result_script = "write_result.pl";

#	Yes, this is a global, for now
$die_on_failed_file_check = config[:die_on_failed_file_check] || false
#$die_on_failed_file_check = config[:die_on_failed_file_check] || true

mailto = config[:mailto] || ''

#	This takes the place of the filenames and paired_end functionality
files = config[:files] || {}


# Possible Errors

raise "File format can either be fastq or fasta" unless( 
	['fasta','fastq'].include?(file_format) )


outdir = "#{Time.now.strftime("%Y%m%d%H%M%S")}.outdir"
FileUtils.mkdir outdir
FileUtils.chdir outdir




# Steps from here


puts "step 1 change fastq files to fasta files"

if( file_format == "fastq")
	files.each_pair do |k,v|
		do_this( "#{fastq2fasta_script} #{v} #{k}lane.fa" );
		file_check( "#{k}lane.fa" );
	end
#if ($file_format eq "fasta") {
else
	if( link_sample_fa_files )
		puts "already fasta format, linking fasta files instead"
		files.each_pair do |k,v|
			do_this( "ln -s #{v} #{k}lane.fa" );
			file_check( "#{k}lane.fa" );
		end
	else
		puts "already fasta format, copy fasta files instead"
		files.each_pair do |k,v|
			do_this( "cp #{v} #{k}lane.fa" );
			file_check( "#{k}lane.fa" );
		end
	end
end



#
#	TODO for some reason, blat doesn't work on the chopped????
#
#	blat34(2007) works, blat35(2012) does not???  Substantial diffs somewhere.
#
puts "step 2 chop reads"
files.each_pair do |k,v|
	do_this( "#{chopreads_script} #{k}lane.fa chopped_#{k}lane.fa #{chop_read_length}" );
	file_check( "chopped_#{k}lane.fa" );
end

puts "step 3 blat chopped reads"
files.each_pair do |k,v|
	do_this( "#{blat_bin} #{blat_reference} -minIdentity=#{minIdentity} chopped_#{k}lane.fa chopped_#{k}lane.psl" );
	file_check( "chopped_#{k}lane.psl", 427 );
end

puts "step 4 find blat out candidate reads"
cmd = "#{blat_out_candidate_script} "
#	files is a hash and the keys are not guaranteed to be sorted
#	sort alphabetically and left is first, right is last (conveniently)
files.keys.sort.each{|k| cmd << "chopped_#{k}lane.psl " }
files.keys.sort.each{|k| cmd << "#{k}lane.fa " }
do_this( cmd );
files.each_pair { |k,v| file_check( "blat_out_candidate_#{k}lane.fa" ) }


puts "step 5 compress raw reads"
files.each_pair do |k,v|
	do_this( "#{compress_script} blat_out_candidate_#{k}lane.fa #{compress_ratio_thrd} > compress_#{k}lane.names" )
	file_check( "compress_#{k}lane.names" )
end


puts "step 6 pull reads from blat_out_candidate fasta files"
cmd = "#{pull_reads_fasta_script} "
#	files is a hash and the keys are not guaranteed to be sorted
#	sort alphabetically and left is first, right is last (conveniently)
files.keys.sort.each{|k| cmd << "compress_#{k}lane.names " }
files.keys.sort.each{|k| cmd << "blat_out_candidate_#{k}lane.fa " }
files.keys.sort.each{|k| cmd << "compress_#{k}lane.fa " }
do_this( cmd );
files.each_pair { |k,v| file_check( "compress_#{k}lane.fa" ) }


puts "step 7 align compressed reads to human genome reference using bowtie"
files.each_pair do |k,v|
	do_this( "#{bowtie_bin} -n #{bowtie_mismatch} -p #{bowtie_threads} -f -S #{bowtie_index_human} compress_#{k}lane.fa compress_#{k}lane.sam" )
	file_check( "compress_#{k}lane.sam" )
	do_this( "#{sam2names_script} compress_#{k}lane.sam bowtie_#{k}lane.names" )
	file_check( "bowtie_#{k}lane.names" )
end
cmd = "#{pull_reads_fasta_script} "
#	files is a hash and the keys are not guaranteed to be sorted
#	sort alphabetically and left is first, right is last (conveniently)
files.keys.sort.each{|k| cmd << "bowtie_#{k}lane.names " }
files.keys.sort.each{|k| cmd << "compress_#{k}lane.fa " }
files.keys.sort.each{|k| cmd << "bowtie_#{k}lane.fa " }
do_this( cmd );
files.each_pair { |k,v| file_check( "bowtie_#{k}lane.fa" ) }

cmd = "#{candidate_non_human_script} "
#	files is a hash and the keys are not guaranteed to be sorted
#	sort alphabetically and left is first, right is last (conveniently)
files.keys.sort.each{|k| cmd << "bowtie_#{k}lane.names " }
do_this( cmd );
file_check( "candidate_non_human.txt" )


nth_iteration = 1;
puts "step 8 #{nth_iteration} iteration"
puts "de novo assembly using Trinity"

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
#	Downgraded to trinityrnaseq_r2011-08-20
#	so can undo these parameter modifications.
#	Hopefully this will change the output as well
#	which would remove the need for the modifications to
#	blastn_cleanup.pl and write_results.pl
#
#	Downgraded Blat from 35 to 34
#	and Blast from 2.2.26 to 2.2.24
#

cmd = "#{trinity_script} --seqType fa --group_pairs_distance #{paired_fragment_length} --min_contig_length #{min_contig_length} --output trinity_output --CPU #{trinity_threads} --bfly_opts \"--stderr\" --JM 1G "
files.each_pair { |k,v| cmd << "--#{k} bowtie_#{k}lane.fa " }
do_this(cmd)
file_check( 'trinity_output/Trinity.fasta' );


do_this( "cp trinity_output/Trinity.fasta Trinity.fasta" );

#
#	This script just joins the sequence on a single line
#	rather that having it in 60 character line segments.
#	I don't know why this is needed, but nevertheless ...
#
do_this( "#{modify_trinity_output_script} Trinity.fasta" );



#	TODO	insert some sort of check


puts "blastn trinity output against human genome"
do_this( "#{blastn_bin} -query=Trinity.fasta -db=#{blastn_index_human} -evalue #{blastn_evalue_thrd} -outfmt 6 > human_contig.txt" )
file_check( 'human_contig.txt' )	#	NOTE  don't know how big an "empty" one is


puts "clean up blastn outputs"
do_this( "#{blastn_cleanup_script} human_contig.txt Trinity.fasta clean_blastn.fa #{similarity_thrd}" )
file_check( 'clean_blastn.fa' )	#	NOTE  don't know how big an "empty" one is



#for (nth_iteration=2; nth_iteration <= iteration; nth_iteration+=1 ) do
for nth_iteration in 2..iteration
	puts "step 8 #{nth_iteration} iteration"
	puts "blat chopped reads"
	files.each_pair do |k,v|
		do_this( "#{blat_bin} clean_blastn.fa -minIdentity=95 #{k}lane.fa #{k}lane.iteration.psl" );
		file_check( "#{k}lane.iteration.psl", 427 );
	end
	puts "find blat out candidate reads"

	cmd = "#{blat_out_candidate_script} "
	files.keys.sort.each{|k| cmd << "#{k}lane.iteration.psl " }
	files.keys.sort.each{|k| cmd << "#{k}lane.fa " }
	do_this( cmd )
	files.each_pair do |k,v|
		file_check( "blat_out_candidate_#{k}lane.fa" )
		do_this( "cp blat_out_candidate_#{k}lane.fa iteration_#{k}lane.fa" );
	end

	print "de novo assembly using Trinity\n";
	do_this( "rm -r trinity_output" );

	cmd = "#{trinity_script} --seqType fa --group_pairs_distance #{paired_fragment_length} --min_contig_length #{min_contig_length} --output trinity_output --CPU #{trinity_threads} --bfly_opts \"--stderr\" --JM 1G "
	files.each_pair { |k,v| cmd << "--#{k} bowtie_#{k}lane.fa " }
	do_this(cmd)
	file_check( 'trinity_output/Trinity.fasta' )

	do_this( "cp trinity_output/Trinity.fasta Trinity.fasta" )
	do_this( "#{modify_trinity_output_script} Trinity.fasta" )

	puts "blastn trinity output against human genome"
	do_this( "#{blastn_bin} -query=Trinity.fasta -db=#{blastn_index_human} -evalue #{blastn_evalue_thrd} -outfmt 6 > human_contig.txt" )
	file_check( 'human_contig.txt' )	#	NOTE  don't know how big an "empty" one is

	puts "clean up blastn outputs"
	do_this( "#{blastn_cleanup_script} human_contig.txt Trinity.fasta clean_blastn.fa #{similarity_thrd}" )
	file_check( 'clean_blastn.fa' );	#	NOTE  don't know how big an "empty" one is

	do_this( "rm -r trinity_output" );
end


puts "step 9 detect species of non human sequences"
puts "blastn trinity output against non-human genome"
do_this( "cp clean_blastn.fa non_human_contig.fa" );
do_this( "#{blastn_bin} -query=non_human_contig.fa -db=#{blastn_index_non_human} -evalue #{blastn_evalue_thrd} -outfmt 6 > non_human_contig_blastn.txt" )
file_check( 'non_human_contig_blastn.txt' );	#	NOTE  don't know how big an "empty" one is


files.each_pair do |k,v|
	do_this( "#{blat_bin} non_human_contig.fa -minIdentity=98 iteration_#{k}lane.fa #{k}lane.psl" );
	file_check( "#{k}lane.psl", 427 )
end

puts "write results"
cmd = "#{write_result_script} non_human_contig.fa non_human_contig_blastn.txt "
files.each_pair { |k,v| cmd << "#{k}lane.psl " }
cmd << output_filename
do_this(cmd)

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

