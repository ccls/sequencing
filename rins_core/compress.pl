#!/usr/bin/perl
#
# This script is used to compress reads used the LZW compression ratio as a hard cutoff
#
# Use: input file, direct to output with > 

use warnings;
use strict;

my $input = $ARGV[0];
my $bits = 8;
my $compress_ratio_thrd = $ARGV[1];

open (IN, "<$input");

my (%dictionary, @dictionary,$namesfastq,$counter);

while (my $line = <IN>){
  chomp $line;
  $counter++;

  $namesfastq = $line;

  $line = <IN>;
  chomp $line;

  my @letters = split (//, $line);
  my $count;
  $dictionary{"a"} = 1;
  $dictionary{"t"} = 1;
  $dictionary{"g"} = 1;
  $dictionary{"c"} = 1;
  push (@dictionary, "a");
  push (@dictionary, "t");
  push (@dictionary, "c");
  push (@dictionary, "g");

  my $end = 1;
  my $startpoint = 0;
  for ($count = 1; $count < (10*($#letters)); $count++){
    my ($start, @tojoin);
    for ($start = $startpoint; $start < $end; $start++){
      push (@tojoin, $letters[$start]);
    }
    my $code = join ("", @tojoin);
    if (exists $dictionary{$code}){
      $end++;
      if ($end == $#letters + 1){
	$dictionary{$bits} = $code;
	push (@dictionary, $bits);
	last;
      }
      next;
    }
    else{
      if ($end == ($#letters + 1)){
	last;
      }
      $dictionary{$code} = $code;
      push (@dictionary, $code);
      $startpoint = ($end-1);
      $end = ($startpoint + 1);
      next;
    }
  }
  my $originalsize = ($#letters + 1)*$bits;
  my $compressed = ($#dictionary + 1)*$bits;
  my $compression_ratio = $compressed/$originalsize;
  if ($compression_ratio > $compress_ratio_thrd){

    $namesfastq =~s/\/[0-9]+$//;
    $namesfastq =~s/^\>//;

    print "$namesfastq\n";
  }
  @dictionary = ();
  %dictionary = ();

}


close (IN);
