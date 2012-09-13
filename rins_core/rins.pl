#!/usr/bin/env perl -w
##!/usr/bin/perl -w

# rins.pl is a tool to Rapidly Identify Nonhuman Sequences
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
push @usage, "Command Line Example: perl rins.pl -c config.txt -o output.txt\n";
push @usage, "  -h, --help        Displays this information\n";
push @usage, "  -o, --output      Output Filename\n";
push @usage, "  -c, --config      Config Filename\n";


my $help;
my $output_filename;
my $config_filename;

system "date";

GetOptions 
(
 'help'        => \$help,
 'output=s'    => \$output_filename,
 'config=s'    => \$config_filename,
);

not defined $help or die @usage;
defined $config_filename or die @usage;
defined $output_filename or die @usage;



my $config = new();
$config->read($config_filename);


# Config values of files and files' configurations

my $file_format = $config->get_value("file_format");
my $pair_end = $config->get_value("pair_end");

my $leftlane_filename  = $config->get_value("leftlane_filename");
my $rightlane_filename = $config->get_value("rightlane_filename");
my $singlelane_filename = $config->get_value("singlelane_filename");

my $raw_read_length = $config->get_value("raw_read_length");
my $chop_read_length = $config->get_value("chop_read_length");
my $minIdentity = $config->get_value("minIdentity");


my $blat_reference = $config->get_value("blat_reference");

my $compress_ratio_thrd = $config->get_value("compress_ratio_thrd");
my $iteration = $config->get_value("iteration");

# executables' configurations

my $blat_bin = $config->get_value("blat_bin");
my $bowtie_bin = $config->get_value("bowtie_bin");
my $bowtie_build_bin = $config->get_value("bowtie_build_bin");
my $bowtie_index_human = $config->get_value("bowtie_index_human");
my $bowtie_threads = $config->get_value("bowtie_threads");
my $bowtie_mismatch = $config->get_value("bowtie_mismatch");

my $trinity_script = $config->get_value("trinity_script");
my $paired_fragment_length = $config->get_value("paired_fragment_length");
my $min_contig_length = $config->get_value("min_contig_length");
my $trinity_threads = $config->get_value("trinity_threads");

my $blastn_bin = $config->get_value("blastn_bin");
my $blastn_index_human = $config->get_value("blastn_index_human");
my $blastn_index_non_human = $config->get_value("blastn_index_non_human");
my $blastn_evalue_thrd = $config->get_value("blastn_evalue_thrd");
my $similarity_thrd = $config->get_value("similarity_thrd");


# scripts' configurations

my $scripts_directory = $config->get_value("scripts_directory");
my $fastq2fasta_script = "$scripts_directory/fastq2fasta.pl";
my $chopreads_script = "$scripts_directory/chopreads.pl";
my $blat_out_candidate_script = "$scripts_directory/blatoutcandidate.pl";
my $compress_script = "$scripts_directory/compress.pl";
my $pull_reads_fasta_script = "$scripts_directory/pull_reads_fasta.pl";
my $sam2names_script = "$scripts_directory/sam2names.pl";
my $modify_trinity_output_script = "$scripts_directory/modify_trinity_output.pl";
my $blastn_cleanup_script = "$scripts_directory/blastn_cleanup.pl";
my $candidate_non_human_script = "$scripts_directory/candidate_non_human.pl";
my $write_result_script = "$scripts_directory/write_result.pl";


my $mailto;
if ($config->has_value("mailto")) {
  $mailto = $config->get_value("mailto");
}


# Possible Errors

if (($file_format ne "fastq") and ($file_format ne "fasta")) {
  die "File format can either be fastq or fasta\n";
}



# Steps from here


print "step 1 change fastq files to fasta files\n";

if ($file_format eq "fastq") {  
  if ($pair_end) {
    system "$fastq2fasta_script $leftlane_filename leftlane.fa";
    system "$fastq2fasta_script $rightlane_filename rightlane.fa";
  }
  else {
    system "$fastq2fasta_script $singlelane_filename singlelane.fa";
  }
}

if ($file_format eq "fasta") {
  print "already fasta format, copy fasta files instead\n";

  if ($pair_end) {
    system "cp $leftlane_filename leftlane.fa";
    system "cp $rightlane_filename rightlane.fa";
  }
  else {
    system "cp $singlelane_filename singlelane.fa";
  }
}




print "step 2 chop reads\n";

if ($pair_end) {
  system "$chopreads_script leftlane.fa chopped_leftlane.fa $chop_read_length";
  system "$chopreads_script rightlane.fa chopped_rightlane.fa $chop_read_length";
}
else {
  system "$chopreads_script singlelane.fa chopped_singlelane.fa $chop_read_length";
}



print "step 3 blat chopped reads\n";

if ($pair_end) {
  system "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_leftlane.fa chopped_leftlane.psl";
  system "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_rightlane.fa chopped_rightlane.psl";
}
else {
  system "$blat_bin $blat_reference -minIdentity=$minIdentity chopped_singlelane.fa chopped_singlelane.psl";
}



print "step 4 find blat out candidate reads\n";

if ($pair_end) {
  system "$blat_out_candidate_script chopped_leftlane.psl chopped_rightlane.psl leftlane.fa rightlane.fa";
}

else {
  system "$blat_out_candidate_script chopped_singlelane.psl singlelane.fa";
}



print "step 5 compress raw reads\n";

if ($pair_end) {
  system "$compress_script blat_out_candidate_leftlane.fa $compress_ratio_thrd > compress_leftlane.names";
  system "$compress_script blat_out_candidate_leftlane.fa $compress_ratio_thrd > compress_rightlane.names";
}

else {
  system "$compress_script blat_out_candidate_singlelane.fa $compress_ratio_thrd > compress_singlelane.names";
}



print "step 6 pull reads from blat_out_candidate fasta files\n";
if ($pair_end) {
  system "$pull_reads_fasta_script compress_leftlane.names compress_rightlane.names blat_out_candidate_leftlane.fa blat_out_candidate_rightlane.fa compress_leftlane.fa compress_rightlane.fa";
}

else {
  system "$pull_reads_fasta_script compress_singlelane.names blat_out_candidate_singlelane.fa compress_singlelane.fa";
}


print "step 7 align compressed reads to human genome reference using bowtie\n";

if ($pair_end) {
  system "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_leftlane.fa compress_leftlane.sam";
  system "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_rightlane.fa compress_rightlane.sam";

  system "$sam2names_script compress_leftlane.sam bowtie_leftlane.names";
  system "$sam2names_script compress_rightlane.sam bowtie_rightlane.names";
  
  system "$pull_reads_fasta_script bowtie_leftlane.names bowtie_rightlane.names compress_leftlane.fa compress_rightlane.fa bowtie_leftlane.fa bowtie_rightlane.fa";
}

else {
  system "$bowtie_bin -n $bowtie_mismatch -p $bowtie_threads -f -S $bowtie_index_human compress_singlelane.fa compress_singlelane.sam";
  system "$sam2names_script compress_singlelane.sam bowtie_singlelane.names";
  system "$pull_reads_fasta_script bowtie_singlelane.names compress_singlelane.fa bowtie_singlelane.fa";
}

if ($pair_end) {
  system "$candidate_non_human_script bowtie_leftlane.names bowtie_rightlane.names";
}
else {
  system "$candidate_non_human_script bowtie_singlelane.names";
}



my $nth_iteration = 1;

print "step 8 $nth_iteration iteration\n";

print "de novo assembly using Trinity\n";

if ($pair_end) {
  system "$trinity_script --seqType fa --left bowtie_leftlane.fa --right bowtie_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"";
}

else {
  system "$trinity_script --seqType fa --single bowtie_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"";
}


system "cp trinity_output/Trinity.fasta Trinity.fasta";
system "$modify_trinity_output_script Trinity.fasta";

print "blastn trinity output against human genome\n";
system "$blastn_bin -query=Trinity.fasta -db=$blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt";


print "clean up blastn outputs\n";
system "$blastn_cleanup_script human_contig.txt Trinity.fasta clean_blastn.fa $similarity_thrd";



for ($nth_iteration=2; $nth_iteration<=$iteration; $nth_iteration++) {

  print "step 8 $nth_iteration iteration\n";  
  print "blat chopped reads\n";

  if ($pair_end) {
    system "$blat_bin clean_blastn.fa -minIdentity=95 leftlane.fa leftlane.iteration.psl";
    system "$blat_bin clean_blastn.fa -minIdentity=95 rightlane.fa rightlane.iteration.psl";
  }
  else {
    system "$blat_bin clean_blastn.fa -minIdentity=95 singlelane.fa singlelane.iteration.psl";
  }


  print "find blat out candidate reads\n";
  
  if ($pair_end) {
    system "$blat_out_candidate_script leftlane.iteration.psl rightlane.iteration.psl leftlane.fa rightlane.fa";

    system "cp blat_out_candidate_leftlane.fa iteration_leftlane.fa";
    system "cp blat_out_candidate_rightlane.fa iteration_rightlane.fa";
    
  }
  else {
    system "$blat_out_candidate_script singlelane.iteration.psl singlelane.fa";
    system "cp blat_out_candidate_singlelane.fa iteration_singlelane.fa";
  }


  print "de novo assembly using Trinity\n";
  system "rm -r trinity_output";

  if ($pair_end) {
    system "$trinity_script --seqType fa --left iteration_leftlane.fa --right iteration_rightlane.fa --paired_fragment_length $paired_fragment_length --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"";
  }  
  else {
    system "$trinity_script --seqType fa --single iteration_singlelane.fa --min_contig_length $min_contig_length --run_butterfly --output trinity_output --CPU $trinity_threads --bfly_opts \"--compatible_path_extension --stderr\"";
  }


  system "cp trinity_output/Trinity.fasta Trinity.fasta";
  system "$modify_trinity_output_script Trinity.fasta";

  print "blastn trinity output against human genome\n";
  system "$blastn_bin -query=Trinity.fasta -db=$blastn_index_human -evalue $blastn_evalue_thrd -outfmt 6 > human_contig.txt";

  print "clean up blastn outputs\n";  
  system "$blastn_cleanup_script human_contig.txt Trinity.fasta clean_blastn.fa $similarity_thrd";
  system "rm -r trinity_output";

}



print "step 9 detect species of non human sequences\n";
print "blastn trinity output against non-human genome\n";
system "cp clean_blastn.fa non_human_contig.fa";
system "$blastn_bin -query=non_human_contig.fa -db=$blastn_index_non_human -evalue $blastn_evalue_thrd -outfmt 6 > non_human_contig_blastn.txt";


if ($pair_end) {
  system "$blat_bin non_human_contig.fa -minIdentity=98 iteration_leftlane.fa leftlane.psl";
  system "$blat_bin non_human_contig.fa -minIdentity=98 iteration_rightlane.fa rightlane.psl";
}
else {
  system "$blat_bin non_human_contig.fa -minIdentity=98 iteration_singlelane.fa singlelane.psl";
}


print "write results\n";

if ($pair_end) {
  system "$write_result_script non_human_contig.fa non_human_contig_blastn.txt leftlane.psl rightlane.psl $output_filename";
}
else {
  system "$write_result_script non_human_contig.fa non_human_contig_blastn.txt singlelane.psl $output_filename";
}

&mailme();

system "date";


sub mailme {
  return if not defined $mailto;
  my $text = "RINS analysis completed";
  print "Attempting to mail $mailto the complete notice\n";
  system "echo '$text' | mail -s 'RINS complete notice' $mailto";
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

sub get_list {
  my $self = shift;
  my $key = shift;
  
  my $index = 1;
  my $key_index = $key.$index;
  my @values;
  while ($self->has_value($key_index))
    {
      push @values, $self->get_value($key_index);
      $index++;
      $key_index = $key.$index;
    }
  
  return @values;
}

sub get_hash {
  my $self = shift;
  my $key = shift;
  
  my $config_values = $self->{"config_values"};
  my $config_filename = $self->{"config_filename"};
  
  defined $config_values and defined $config_filename or die "Error: config not read\n";
  
  if (not defined $config_values->{$key})
    {
      die "Error: no value for $key in config file $config_filename\n";
    }
  
  my %values;
  foreach my $value (split /,/, $config_values->{$key})
    {
      $value =~ s/\s//g;
      $values{$value} = 1;
    }
  
  return %values;
}



