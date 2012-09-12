#!/usr/bin/perl
#
# pull out the names of non-mappable sequence
#
#############################################


my $samfile = $ARGV[0];
my $namefile = $ARGV[1];
my $line;
my @data;

open (in, "<$samfile");
open (out, ">$namefile");

while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  if ($data[2] eq "\*") {
    my $name = $data[0];
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;
    print out "$name\n";
  }
}

close in;
close out;
