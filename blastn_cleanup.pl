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

#	blastn -outfmt 6 ...
#	'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
#	comp430_c0_seq1	chrX	97.24	398	10	1	1	398	108184401	108184005	0.0	 673comp430_c0_seq1
  my @data = split /\t/;
  my $numbermatch = $data[3];	#	never seems to be used?
  my $evalue = $data[10];			#	never used either?
  my $bit_score = $data[-1];	#	$data[11] same thing


#	It is my opinion that the new lines are very different from what
#	was originally anticipated.
#	$numbermatch is set to data[3], what I am using as length, 
#		but then it is never used????
#
#	I have very little faith in this


#print "$data[0]\n";
#  my ($front, $back) = split (/len[:]/, $data[0]);
#print "$front .... $back\n";	#	there is no "len[:]" in $data[0] !!!
#  my @length = split (/[_]/, $back);
#  my $percent = $bit_score/$length[0];
  my $percent = $bit_score/$data[3];	#	what if data[3] is actually 0?
#	$percent is not actually a percent.  just a ratio.
#	similarity_thrd is usually 0.8
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
