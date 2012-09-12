#!/usr/bin/perl


my $input = $ARGV[0];
my $fasta_input = $ARGV[1];
my $output = $ARGV[2];
my $similarity_thrd = $ARGV[3];

my %flag = ();
my $line;

open (IN, "<$input");
while (<IN>){
  chomp;
  my @data = split /\t/;
  my $numbermatch = $data[3];
  my $evalue = $data[10];
  my $bit_score = $data[-1];

  my ($front, $back) = split (/len[:]/, $data[0]);
  my @length = split (/[_]/, $back);
  my $percent = $bit_score/$length[0];

  if ($percent >= $similarity_thrd){
    $flag{$data[0]} = 1;
  }
}

close IN;

open (in, "<$fasta_input");
open (out, ">$output");
while ($line=<in>) {
  chomp $line;

  my $name = $line;
  $name =~s/\/[0-9]+$//;
  $name =~s/^\>//;
  
  if (!($flag{$name})) {
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
