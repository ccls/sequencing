#!/usr/bin/env perl
##!/usr/bin/perl



# blastn_cleanup_script human_contig.txt Trinity.fasta clean_blastn.fa 0.8


my $input = $ARGV[0];            #	human_contig.txt
my $fasta_input = $ARGV[1];      #	Trinity.fasta
my $output = $ARGV[2];           #	clean_blastn.fa
my $similarity_thrd = $ARGV[3];  #	0.8

my %flag = ();
my $line;

open (IN, "<$input");
while (<IN>){
  chomp;

#
#	I think that the format of human_contig.txt has changed.
#	This is the format of a given output file, 
#
#comp7_c0_seq1_FPKM_all:46727.836_FPKM_rel:46727.836_len:377_path:[0,586,641,781,805,811]	chr3	99.26	136	0	1	242	377	197266575	197266441	2e-62	 244
#comp7_c0_seq1_FPKM_all:46727.836_FPKM_rel:46727.836_len:377_path:[0,586,641,781,805,811]	chr3	100.00	83	0	0	167	249	197269633	197269551	4e-35	 154


#	blastn -outfmt 6 ...
#	'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
#	comp430_c0_seq1	chrX	97.24	398	10	1	1	398	108184401	108184005	0.0	 673comp430_c0_seq1
  my @data = split /\t/;
#  my $numbermatch = $data[3];	#	never seems to be used?
#  my $evalue = $data[10];			#	never used either?
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

	my $length = $data[3];
  my $percent = $bit_score/$length;	#	what if length is actually 0?
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
