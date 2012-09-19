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

config_filename = 'config.yml'
output_filename = 'results.txt'
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
pair_end    = config[:pair_end]
link_sample_fa_files = config[:link_sample_fa_files] || false
leftlane_filename	  = config[:leftlane_filename]
rightlane_filename  = config[:rightlane_filename]
singlelane_filename = config[:singlelane_filename]
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

die_on_failed_file_check = config[:die_on_failed_file_check] || false
#die_on_failed_file_check = config[:die_on_failed_file_check] || true

mailto = config[:mailto] || ''



# Possible Errors

raise "File format can either be fastq or fasta" unless( 
	['fasta','fastq'].include?(file_format) )

__END__




#	make and use an outdir based on time now.
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;
my $outdir = sprintf( "%4d%02d%02d%02d%02d%02d.outdir",
	$year+1900, $mon+1, $mday, $hour, $min, $sec );
mkdir $outdir;
chdir $outdir;





__END__


# Steps from here


print "step 1 change fastq files to fasta files\n";
if ($file_format eq "fastq") {	
	if ($pair_end) {
		do_this( "$fastq2fasta_script $leftlane_filename leftlane.fa" );
		do_this( "$fastq2fasta_script $rightlane_filename rightlane.fa" );
		file_check( 'leftlane.fa' );
		file_check( 'rightlane.fa' );
	} else {
		do_this( "$fastq2fasta_script $singlelane_filename singlelane.fa" );
		file_check( 'singlelane.fa' );
	}
}

if ($file_format eq "fasta") {
	if( $link_sample_fa_files =~ /t/i ){
		print "already fasta format, linking fasta files instead\n";
		if ($pair_end) {
			do_this( "ln -s $leftlane_filename leftlane.fa" );
			do_this( "ln -s $rightlane_filename rightlane.fa" );
			file_check( 'leftlane.fa' );
			file_check( 'rightlane.fa' );
		} else {
			do_this( "ln -s $singlelane_filename singlelane.fa" );
			file_check( 'singlelane.fa' );
		}
	} else {
		print "already fasta format, copy fasta files instead\n";
		if ($pair_end) {
			do_this( "cp $leftlane_filename leftlane.fa" );
			do_this( "cp $rightlane_filename rightlane.fa" );
			file_check( 'leftlane.fa' );
			file_check( 'rightlane.fa' );
		} else {
			do_this( "cp $singlelane_filename singlelane.fa" );
			file_check( 'singlelane.fa' );
		}
	}
}

print "step 2 chop reads\n";

#
#	TODO for some reason, blat doesn't work on the chopped????
#
#	blat34(2007) works, blat35(2012) does not???  Substantial diffs somewhere.
#
if ($pair_end) {
	do_this( "$chopreads_script leftlane.fa chopped_leftlane.fa $chop_read_length" );
#	do_this( "cp $leftlane_filename chopped_leftlane.fa" );
	do_this( "$chopreads_script rightlane.fa chopped_rightlane.fa $chop_read_length" );
#	do_this( "cp $rightlane_filename chopped_rightlane.fa" );
	file_check( 'chopped_leftlane.fa' );
	file_check( 'chopped_rightlane.fa' );

} else {
	do_this( "$chopreads_script singlelane.fa chopped_singlelane.fa $chop_read_length" );
#	do_this( "cp $singlelane_filename chopped_singlelane.fa" );
	file_check( 'chopped_singlelane.fa' );
}


print "step 3 blat chopped reads\n";

if ($pair_end) {
	do_this( "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_leftlane.fa chopped_leftlane.psl" );
	do_this( "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_rightlane.fa chopped_rightlane.psl" );
	file_check( 'chopped_leftlane.psl', 427 );
	file_check( 'chopped_rightlane.psl', 427 );
} else {
	do_this( "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_singlelane.fa chopped_singlelane.psl" );
	file_check( 'chopped_singlelane.psl', 427 );
}

#
#
#	These psl files are empty when generated from the chopped_*lane.fa files
#
#


print "step 4 find blat out candidate reads\n";

if ($pair_end) {
	do_this( "$blat_out_candidate_script chopped_leftlane.psl chopped_rightlane.psl leftlane.fa rightlane.fa" );
	file_check( 'blat_out_candidate_leftlane.fa' );
	file_check( 'blat_out_candidate_rightlane.fa' );
} else {
	do_this( "$blat_out_candidate_script chopped_singlelane.psl singlelane.fa" );
	file_check( 'blat_out_candidate_singlelane.fa' );
}

print "step 5 compress raw reads\n";

if ($pair_end) {
	do_this( "$compress_script blat_out_candidate_leftlane.fa $compress_ratio_thrd > compress_leftlane.names" );
	file_check( 'compress_leftlane.names' );
	do_this( "$compress_script blat_out_candidate_rightlane.fa $compress_ratio_thrd > compress_rightlane.names" );
	file_check( 'compress_rightlane.names' );
} else {
	do_this( "$compress_script blat_out_candidate_singlelane.fa $compress_ratio_thrd > compress_singlelane.names" );
	file_check( 'compress_singlelane.names' );
}


print "step 6 pull reads from blat_out_candidate fasta files\n";
if ($pair_end) {
	do_this( "$pull_reads_fasta_script compress_leftlane.names compress_rightlane.names blat_out_candidate_leftlane.fa blat_out_candidate_rightlane.fa compress_leftlane.fa compress_rightlane.fa" );
	file_check( 'compress_leftlane.fa' );
	file_check( 'compress_rightlane.fa' );
} else {
	do_this( "$pull_reads_fasta_script compress_singlelane.names blat_out_candidate_singlelane.fa compress_singlelane.fa" );
	file_check( 'compress_singlelane.fa' );
}


print "step 7 align compressed reads to human genome reference using bowtie\n";

if ($pair_end) {
	do_this( "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_leftlane.fa compress_leftlane.sam" );
	do_this( "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_rightlane.fa compress_rightlane.sam" );
	file_check( 'compress_leftlane.sam' );
	file_check( 'compress_rightlane.sam' );

	do_this( "$sam2names_script compress_leftlane.sam bowtie_leftlane.names" );
	do_this( "$sam2names_script compress_rightlane.sam bowtie_rightlane.names" );
	file_check( 'bowtie_leftlane.names' );
	file_check( 'bowtie_rightlane.names' );
	
	do_this( "$pull_reads_fasta_script bowtie_leftlane.names bowtie_rightlane.names compress_leftlane.fa compress_rightlane.fa bowtie_leftlane.fa bowtie_rightlane.fa" );
	file_check( 'bowtie_leftlane.fa' );
	file_check( 'bowtie_rightlane.fa' );
} else {
	do_this( "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_singlelane.fa compress_singlelane.sam" );
	file_check( 'compress_singlelane.sam' );
	do_this( "$sam2names_script compress_singlelane.sam bowtie_singlelane.names" );
	file_check( 'bowtie_singlelane.names' );
	do_this( "$pull_reads_fasta_script bowtie_singlelane.names compress_singlelane.fa bowtie_singlelane.fa" );
	file_check( 'bowtie_singlelane.fa' );
}

if ($pair_end) {
	do_this( "$candidate_non_human_script bowtie_leftlane.names bowtie_rightlane.names" );
}
else {
	do_this( "$candidate_non_human_script bowtie_singlelane.names" );
}
file_check( 'candidate_non_human.txt' );



my $nth_iteration = 1;

print "step 8 $nth_iteration iteration\n";

print "de novo assembly using Trinity\n";



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



if ($pair_end) {
	#	trinityrnaseq_r2011-08-20
	#	do_this( "$trinity_script --seqType fa --left bowtie_leftlane.fa --right bowtie_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\"" );

	#	original
	#	do_this( "$trinity_script --seqType fa --left bowtie_leftlane.fa --right bowtie_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"" );

	#	latest
	do_this( "$trinity_script --seqType fa --left bowtie_leftlane.fa --right bowtie_rightlane.fa --group_pairs_distance $paired_fragment_length --min_contig_length $min_contig_length --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\" --JM 1G" );
} else {
	#	trinityrnaseq_r2011-08-20
	#	do_this( "$trinity_script --seqType fa --single bowtie_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\"" );

	#	original
	#	do_this( "$trinity_script --seqType fa --single bowtie_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"" );

	#	latest
	do_this( "$trinity_script --seqType fa --single bowtie_singlelane.fa --min_contig_length $min_contig_length --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\" --JM 1G" );
}

file_check( 'trinity_output/Trinity.fasta' );

do_this( "cp trinity_output/Trinity.fasta Trinity.fasta" );

#
#	This script just joins the sequence on a single line
#	rather that having it in 60 character line segments.
#	I don't know why this is needed, but nevertheless ...
#
do_this( "$modify_trinity_output_script Trinity.fasta" );



#	TODO	insert some sort of check


print "blastn trinity output against human genome\n";
do_this( "$blastn_bin -query=Trinity.fasta -db=$blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt" );

file_check( 'human_contig.txt' );	#	NOTE  don't know how big an "empty" one is

print "clean up blastn outputs\n";
do_this( "$blastn_cleanup_script human_contig.txt Trinity.fasta clean_blastn.fa $similarity_thrd" );

file_check( 'clean_blastn.fa' );	#	NOTE  don't know how big an "empty" one is

for ($nth_iteration=2; $nth_iteration<=$iteration; $nth_iteration++) {

	print "step 8 $nth_iteration iteration\n";	
	print "blat chopped reads\n";

	if ($pair_end) {
		do_this( "$blat_bin clean_blastn.fa -minIdentity=95 leftlane.fa leftlane.iteration.psl" );
		do_this( "$blat_bin clean_blastn.fa -minIdentity=95 rightlane.fa rightlane.iteration.psl" );
		file_check( 'leftlane.iteration.psl', 427 );
		file_check( 'rightlane.iteration.psl', 427 );
	}
	else {
		do_this( "$blat_bin clean_blastn.fa -minIdentity=95 singlelane.fa singlelane.iteration.psl" );
		file_check( 'singlelane.iteration.psl', 427 );
	}

	print "find blat out candidate reads\n";
	
	if ($pair_end) {
		do_this( "$blat_out_candidate_script leftlane.iteration.psl rightlane.iteration.psl leftlane.fa rightlane.fa" );
		file_check( 'blat_out_candidate_leftlane.fa' );
		file_check( 'blat_out_candidate_rightlane.fa' );

		do_this( "cp blat_out_candidate_leftlane.fa iteration_leftlane.fa" );
		do_this( "cp blat_out_candidate_rightlane.fa iteration_rightlane.fa" );
		
	}
	else {
		do_this( "$blat_out_candidate_script singlelane.iteration.psl singlelane.fa" );
		file_check( 'blat_out_candidate_singlelane.fa' );
		do_this( "cp blat_out_candidate_singlelane.fa iteration_singlelane.fa" );
	}

	print "de novo assembly using Trinity\n";
	do_this( "rm -r trinity_output" );

	if ($pair_end) {
		#	trinityrnaseq_r2011-08-20
		#	do_this( "$trinity_script --seqType fa --left iteration_leftlane.fa --right iteration_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\"" );

		#	original
		#	do_this( "$trinity_script --seqType fa --left iteration_leftlane.fa --right iteration_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"" );

		#	latest
		do_this( "$trinity_script --seqType fa --left iteration_leftlane.fa --right iteration_rightlane.fa --group_pairs_distance $paired_fragment_length --min_contig_length $min_contig_length --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\" --JM 1G" );
	}	
	else {
		#	trinityrnaseq_r2011-08-20
		#	do_this( "$trinity_script --seqType fa --single iteration_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\"" );

		#	original
		#	do_this( "$trinity_script --seqType fa --single iteration_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"" );

		#	latest
		do_this( "$trinity_script --seqType fa --single iteration_singlelane.fa --min_contig_length $min_contig_length --output trinity_output --CPU $trinity_threads --bfly_opts \"--stderr\" --JM 1G" );
	}

	file_check( 'trinity_output/Trinity.fasta' );

	do_this( "cp trinity_output/Trinity.fasta Trinity.fasta" );
	do_this( "$modify_trinity_output_script Trinity.fasta" );

	print "blastn trinity output against human genome\n";
	do_this( "$blastn_bin -query=Trinity.fasta -db=$blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt" );
	file_check( 'human_contig.txt' );	#	NOTE  don't know how big an "empty" one is

	print "clean up blastn outputs\n";	
	do_this( "$blastn_cleanup_script human_contig.txt Trinity.fasta clean_blastn.fa $similarity_thrd" );
	file_check( 'clean_blastn.fa' );	#	NOTE  don't know how big an "empty" one is

	do_this( "rm -r trinity_output" );
}





print "step 9 detect species of non human sequences\n";
print "blastn trinity output against non-human genome\n";
do_this( "cp clean_blastn.fa non_human_contig.fa" );
do_this( "$blastn_bin -query=non_human_contig.fa -db=$blastn_index_non_human -evalue $blastn_evalue_thrd -outfmt 6 > non_human_contig_blastn.txt" );

file_check( 'non_human_contig_blastn.txt' );	#	NOTE  don't know how big an "empty" one is

if ($pair_end) {
	do_this( "$blat_bin non_human_contig.fa -minIdentity=98 iteration_leftlane.fa leftlane.psl" );
	do_this( "$blat_bin non_human_contig.fa -minIdentity=98 iteration_rightlane.fa rightlane.psl" );
	file_check( 'leftlane.psl', 427 );
	file_check( 'rightlane.psl', 427 );
}
else {
	do_this( "$blat_bin non_human_contig.fa -minIdentity=98 iteration_singlelane.fa singlelane.psl" );
	file_check( 'singlelane.psl', 427 );
}

print "write results\n";

if ($pair_end) {
	do_this( "$write_result_script non_human_contig.fa non_human_contig_blastn.txt leftlane.psl rightlane.psl $output_filename" );
}
else {
	do_this( "$write_result_script non_human_contig.fa non_human_contig_blastn.txt singlelane.psl $output_filename" );
}



#	TODO	insert some sort of check


&mailme();

system "date";

######################################################################


sub file_check {
	my ( $filename, $empty_size ) = @_;
	$empty_size ||= 0;	# usually 0, but the psl files are 427 (header only)

#	STDERR is not logged when using " | tee -a log"

	my $msg = '';
	unless( -e $filename ){
		$msg = "$filename not created";
		if( $die_on_failed_file_check ){
			die $msg;
		} else {
			print "$msg\n";
		}
	}
	#	Need to check existance too because if not dying on failure
	#	could have passed above check.  Then '-s $filename' is '' raising ...
	#	Use of uninitialized value in numeric le (<=) at ...
	if( ( -e $filename ) && ( -s $filename <= $empty_size ) ){
		$msg = "$filename empty ( <= $empty_size )";
		if( $die_on_failed_file_check ){
			die $msg
		} else {
			print "$msg\n";
		}
	}
}

sub do_this {
	my ( $command ) = @_;
	print ( "Executing ...\n" );
	print ( "$command\n" );

#	backticks cache the output so you see nothing until its done.
#	@result would be that output that COULD be printed.
#	my @result = `$command`;	

#	$exit_status here is the same as the built in $?
#	so long as another system command hasn't been run.
#	my $exit_status = system($command);
	system($command);

	print ( "\n$command failed with $?\n\n" ) if ( $? );
} 








sub mailme {
	return if not defined $mailto;
	my $text = "RINS analysis completed";
	print "Attempting to mail $mailto the complete notice\n";
#	system "echo '$text' | mail -s 'RINS complete notice' $mailto";
}

sub new {
	my $self = {};
	bless($self);
	return $self;
}

sub read {
	my $self = shift;
	my $config_filename = shift;
	my %config_values;
	
	open CFG, $config_filename or die "Error: Unable to open $config_filename\n";
	while (<CFG>) {
		chomp;
		/^\s*([^=\s]+)\s*=\s*(.*)$/;
			
		my $key = $1;
		my $value = $2;
			
		next if not defined $key;
		next if not defined $value;
		
		$config_values{$key} = $value;

	}
	close CFG;
	
#
#	I don't think that this is used?  Looks like some sort of self referential
#
#	foreach my $key (keys %config_values) {
#		while ($config_values{$key} =~ /\$\(([^)]+)\)/) {
#			my $other_key = $1;
#		
#			if (not defined $config_values{$other_key}) {
#				die "Error: no value for $other_key in config file $config_filename\n";
#			}
#		
#			$config_values{$key} =~ s/\$\($other_key\)/$config_values{$other_key}/;
#		}
#	}

	$self->{"config_values"} = \%config_values;
	$self->{"config_filename"} = $config_filename;
}

sub has_value {
	my $self = shift;
	my $key = shift;
	
	my $config_values = $self->{"config_values"};
	my $config_filename = $self->{"config_filename"};
	
	defined $config_values and defined $config_filename or die "Error: config not read\n";
	
	return defined $config_values->{$key};
}

sub get_value {
	my $self = shift;
	my $key = shift;
	
	my $config_values = $self->{"config_values"};
	my $config_filename = $self->{"config_filename"};
	
	defined $config_values and defined $config_filename or die "Error: config not read\n";
		
	return $config_values->{$key};
}


#
#	A list in the config file looks like ...
#	mylist0 = apple
#	mylist1 = orange
#	mylist2 = banana
#
#sub get_list {
#	my $self = shift;
#	my $key = shift;
#	
##	originally started with 1, but the array index starts with 0
##	my $index = 1;	
##	avoid the discrepancy and start with 0 too.
##	(never used anyway)
#
#	my $index = 0;
#	my $key_index = $key.$index;
#	my @values;
#	while ($self->has_value($key_index)) {
#		push @values, $self->get_value($key_index);
#		$index++;
#		$key_index = $key.$index;
#	}
#	
#	return @values;
#}

#sub get_hash {
#	my $self = shift;
#	my $key = shift;
#	
#	my $config_values = $self->{"config_values"};
#	my $config_filename = $self->{"config_filename"};
#	
#	defined $config_values and defined $config_filename or die "Error: config not read\n";
#	
#	if (not defined $config_values->{$key})
#		{
#			die "Error: no value for $key in config file $config_filename\n";
#		}
#	
#	my %values;
#	foreach my $value (split /,/, $config_values->{$key})
#		{
#			$value =~ s/\s//g;
#			$values{$value} = 1;
#		}
#	
#	return %values;
#}

