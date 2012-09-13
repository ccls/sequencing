#!/usr/bin/env perl -w
##!/usr/bin/perl -w

# rins.pl is a tool to Rapidly Identify Nonhuman Sequences
#	(This is a CCLS modified version of the Stanford rins.pl by Jakeb)
# 
# script is written by Kun Qu and Aparna Bhaduri in Stanford Dermatology

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use File::Spec;


my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Run the RINS pipeline to identify nonhuman sequences.\n";
push @usage, "Command Line Example: rins.pl -c config.txt -o output.txt\n";
push @usage, "	-h, --help    Displays this information\n";
push @usage, "	-c, --config  Config Filename (default: config.txt)\n";
push @usage, "	-o, --output  Output Filename (default: results.txt)\n";


my $help;
my $output_filename;
my $config_filename;

system "date";

GetOptions 
(
 'help'				=> \$help,
 'output=s'		=> \$output_filename,
 'config=s'		=> \$config_filename,
);

not defined $help or die @usage;
$config_filename ||= 'config.txt';
#defined $config_filename or die @usage;
$output_filename ||= 'results.txt';
#defined $output_filename or die @usage;

my $config = new();
$config->read($config_filename);


# Config values of files and files' configurations

my $file_format = $config->get_value("file_format");
my $pair_end = $config->get_value("pair_end");

my $leftlane_filename	= $config->get_value("leftlane_filename");
my $rightlane_filename = $config->get_value("rightlane_filename");
my $singlelane_filename = $config->get_value("singlelane_filename");

my $raw_read_length = $config->get_value("raw_read_length") || 100;
my $chop_read_length = $config->get_value("chop_read_length") || 25;
my $minIdentity = $config->get_value("minIdentity") || 80;

my $blat_reference = $config->get_value("blat_reference");

my $compress_ratio_thrd = $config->get_value("compress_ratio_thrd") || 0.5;
my $iteration = $config->get_value("iteration") || 2;

# executables' configurations

my $blat_bin = $config->get_value("blat_bin") || 'blat';
my $bowtie_bin = $config->get_value("bowtie_bin") || 'bowtie';
my $bowtie_build_bin = $config->get_value("bowtie_build_bin") || 'bowtie-build';
my $bowtie_index_human = $config->get_value("bowtie_index_human");
my $bowtie_threads = $config->get_value("bowtie_threads") || 6;
my $bowtie_mismatch = $config->get_value("bowtie_mismatch") || 3;

my $trinity_script = $config->get_value("trinity_script") || 'Trinity.pl';
my $paired_fragment_length = $config->get_value("paired_fragment_length") || 300;
my $min_contig_length = $config->get_value("min_contig_length") || 300;
my $trinity_threads = $config->get_value("trinity_threads") || 6;

my $blastn_bin = $config->get_value("blastn_bin") || 'blastn';
my $blastn_index_human = $config->get_value("blastn_index_human");
my $blastn_index_non_human = $config->get_value("blastn_index_non_human");
my $blastn_evalue_thrd = $config->get_value("blastn_evalue_thrd") || 0.05;
my $similarity_thrd = $config->get_value("similarity_thrd") || 0.8;


# scripts' configurations (scripts are all in my path, so dir unnecessary)

#my $scripts_directory = $config->get_value("scripts_directory");
#my $fastq2fasta_script = "$scripts_directory/fastq2fasta.pl";
#my $chopreads_script = "$scripts_directory/chopreads.pl";
#my $blat_out_candidate_script = "$scripts_directory/blatoutcandidate.pl";
#my $compress_script = "$scripts_directory/compress.pl";
#my $pull_reads_fasta_script = "$scripts_directory/pull_reads_fasta.pl";
#my $sam2names_script = "$scripts_directory/sam2names.pl";
#my $modify_trinity_output_script = "$scripts_directory/modify_trinity_output.pl";
#my $blastn_cleanup_script = "$scripts_directory/blastn_cleanup.pl";
#my $candidate_non_human_script = "$scripts_directory/candidate_non_human.pl";
#my $write_result_script = "$scripts_directory/write_result.pl";

my $fastq2fasta_script = "fastq2fasta.pl";
my $chopreads_script = "chopreads.pl";
my $blat_out_candidate_script = "blatoutcandidate.pl";
my $compress_script = "compress.pl";
my $pull_reads_fasta_script = "pull_reads_fasta.pl";
my $sam2names_script = "sam2names.pl";
my $modify_trinity_output_script = "modify_trinity_output.pl";
my $blastn_cleanup_script = "blastn_cleanup.pl";
my $candidate_non_human_script = "candidate_non_human.pl";
my $write_result_script = "write_result.pl";

my $die_on_failed_file_check = $config->get_value("die_on_failed_file_check") || 0;# false
#my $die_on_failed_file_check = $config->get_value("die_on_failed_file_check") || 1;# true


my $mailto;
if ($config->has_value("mailto")) {
	$mailto = $config->get_value("mailto");
}


# Possible Errors

if (($file_format ne "fastq") and ($file_format ne "fasta")) {
	die "File format can either be fastq or fasta\n";
}


#	make and use an outdir based on time now.
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime time;
my $outdir = sprintf( "%4d%02d%02d%02d%02d%02d.outdir",
	$year+1900, $mon+1, $mday, $hour, $min, $sec );
mkdir $outdir;
chdir $outdir;



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
	print "already fasta format, copy fasta files instead\n";
	if ($pair_end) {
		do_this( "cp $leftlane_filename leftlane.fa" );
		do_this( "cp $rightlane_filename rightlane.fa" );
		file_check( 'leftlane.fa' );
		file_check( 'rightlane.fa' );
	}
	else {
		do_this( "cp $singlelane_filename singlelane.fa" );
		file_check( 'singlelane.fa' );
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
	while (<CFG>)
		{
			chomp;
			/^\s*([^=\s]+)\s*=\s*(.*)$/;
			
			my $key = $1;
			my $value = $2;
			
			next if not defined $key;
			next if not defined $value;
			
			$config_values{$key} = $value;
		}
	close CFG;
	
	foreach my $key (keys %config_values)
		{
			while ($config_values{$key} =~ /\$\(([^)]+)\)/)
	{
		my $other_key = $1;
		
		if (not defined $config_values{$other_key})
			{
				die "Error: no value for $other_key in config file $config_filename\n";
			}
		
		$config_values{$key} =~ s/\$\($other_key\)/$config_values{$other_key}/;
	}
		}
	
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

#sub get_list {
#	my $self = shift;
#	my $key = shift;
#	
#	my $index = 1;
#	my $key_index = $key.$index;
#	my @values;
#	while ($self->has_value($key_index))
#		{
#			push @values, $self->get_value($key_index);
#			$index++;
#			$key_index = $key.$index;
#		}
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



