#!/usr/bin/perl

# This script is used to change .fastq file to .fasta file

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

if ($#ARGV ne 1) {
  print "command line example: perl fastq2fasta.pl input.fastq output.fasta\n";
  exit;
}

$infile = $ARGV[0];
$outfile = $ARGV[1];


open (in, "<$infile");
open (out, ">$outfile");
while ($line=<in>) {
  
  $count ++;
  if ($count%1000000 eq 0) {
    print "$count\n";
  }

  chomp $line;
  $line=~s/\s+/\_/g;
  print out ">$line\n";
  $line = <in>;
  print out $line;
  $line = <in>;
  $line = <in>;
}

close in;
close out;
