#!/usr/bin/env perl
##!/usr/bin/perl
#
###########################################################################
#
# this script is used to get unique sequence names from blat out .psl file
#
############################################################################


my %flag;


# pair_end

if ($#ARGV eq 3) {

  for (my $i=0; $i<=1; $i++) {
    open (in, "<$ARGV[$i]");
    my $line = <in>;
    while ($line=<in>) {
      chomp $line;
      my @data = split /\t/, $line;
      my $name = $data[9];
      $name =~s/\/[0-9]+$//;
      $flag{$name} = 1;
    }
    close in;
  }
    
  open (in, "<$ARGV[2]");
  open (out, ">blat_out_candidate_leftlane.fa");

  while ($line=<in>) {
    chomp $line;
    $name = $line;
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;

    if ($flag{$name}) {
      print out "$line\n";
      $line = <in>;
      print out $line;
    }
    else {
      $line = <in>;
    }
  }

  close in;
  close out;

  open (in, "<$ARGV[3]");
  open (out, ">blat_out_candidate_rightlane.fa");

  while ($line=<in>) {
    chomp $line;
    $name = $line;
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;

    if ($flag{$name}) {
      print out "$line\n";
      $line = <in>;
      print out $line;
    }
    else {
      $line = <in>;
    }
  }

  close in;
  close out;

}


# single_end

if ($#ARGV eq 1) {

  open (in, "<$ARGV[0]");
  my $line = <in>;
  while ($line=<in>) {
    chomp $line;
    my @data = split /\t/, $line;
    my $name = $data[9];
    $name =~s/\/[0-9]+$//;
    $flag{$name} = 1;
  }
  close in;
  
    
  open (in, "<$ARGV[1]");
  open (out, ">blat_out_candidate_singlelane.fa");
  
  while ($line=<in>) {
    chomp $line;
    $name = $line;
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;
    
    if ($flag{$name}) {
      print out "$line\n";
      $line = <in>;
      print out $line;
    }
    else {
      $line = <in>;
    }
  }
  
  close in;
  close out;

}

