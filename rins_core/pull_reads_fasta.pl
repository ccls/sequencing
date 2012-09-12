#!/usr/bin/perl
#
###########################################################################
#
# this script is used to get unique sequence names from blat out .psl file
#
############################################################################


my %flag = ();


# pair_end

if ($#ARGV eq 5) {

  for (my $i=0; $i<=1; $i++) {
    open (in, "<$ARGV[$i]");
    while ($line=<in>) {
      chomp $line;
      $name = $line;
      $flag{$name} = 1;
    }
    close in;
  }

  open (in, "<$ARGV[2]");
  open (out, ">$ARGV[4]");

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
  open (out, ">$ARGV[5]");

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

if ($#ARGV eq 2) {

  open (in, "<$ARGV[0]");
  while ($line=<in>) {
    chomp $line;
    $name = $line;
    $flag{$name} = 1;
  }
  close in;
  
    
  open (in, "<$ARGV[1]");
  open (out, ">$ARGV[2]");
  
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

