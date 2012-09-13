#!/usr/bin/env perl
##!/usr/bin/perl
#
# this script is used to modify Trinity.fasta file
#
###################################################


open (in, "<$ARGV[0]");
open (out, ">tmp.txt");

my $line;
my $name;
my %seq;


while ($line=<in>) {
  chomp $line;

  if ($line=~/^\>/) {
    $name = $line;
  }
  else {
    $seq{$name} = $seq{$name}."$line";
  }
}

my @names = keys %seq;
foreach $name (@names) {
  print out "$name\n";
  print out "$seq{$name}\n";
}

close in;
close out;

system "mv tmp.txt $ARGV[0]";
