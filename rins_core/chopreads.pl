#!/usr/bin/perl
#
# This script is used to cut reads in fasta format into 25mers
#
# Use: input fasta file, ouput chopped_fasta files
#
###############################################################

use warnings;
use strict;

my $fasta = $ARGV[0];
my $out = $ARGV[1];
my $chop_read_length = $ARGV[2];


open (FASTA, "<$fasta");
open (OUT, ">$out");

my ($count, %namesfasta, $title);

while (my $line = <FASTA>) {
  chomp $line;
  $count++;
  if ($count%1000000 == 0){
    print "$count\n";
  }
  if ($count%2 == 1) {
    $title = $line;
  }
  
  if ($count%2 == 0) {
    my @sequence = split (//, $line);
    my $totalbases = $#sequence + 1;
    my $number = int($totalbases/$chop_read_length);

    for (my $i=1; $i<=$number; $i++){

      my $start = ($i-1)*$chop_read_length;
      my $end = $i*$chop_read_length - 1;
      my $seq =join ("", @sequence[$start..$end]);
      print OUT "$title\n$seq\n";

    }

    if (($totalbases - $number*$chop_read_length)>=10) {

      my $start = $#sequence - $chop_read_length + 1;
      my $end = $#sequence;
      my $seq = join ("", @sequence[$start..$end]);
      print OUT "$title\n$seq\n";

    }
  }
} 

close (FASTA);
close (OUT);

