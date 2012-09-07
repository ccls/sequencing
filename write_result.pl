#!/usr/bin/perl
#
# used to clean up data and print out results
#
#############################################


my $contig_fa_file = $ARGV[0];
my $contig_blastn_file = $ARGV[1];
my $output_file = $ARGV[-1];
my @psl_files;
for (my $i=2; $i<$#ARGV; $i++) {
  push (@psl_files, $ARGV[$i]);
}


my $line;
my $name;
my $psl_file;
my %seq = ();
my %count = ();


open (in, "<$contig_fa_file");
while ($line=<in>) {
  chomp $line;
  $name = $line;
  $name =~s/\/[0-9]+$//;
  $name =~s/^\>//;

  $line = <in>;
  chomp $line;
  $seq{$name} = $line;

}
close in;


foreach $psl_file (@psl_files) {
  open (in, "<$psl_file");

  for ($i=1;$i<=5;$i++) {
    $line = <in>;
  }

  while ($line=<in>) {
    chomp $line;
    my @data = split /\t/, $line;
    $count{$data[13]} ++; 
  }
  close in;
}


open (in, "<$contig_blastn_file");
open (out, ">$output_file");

print out "contig_name\tnumber_of_raw_reads_fall_on_this_contig\tnon_human_species\tE-value\tbit_score\tcontig_sequence\n";

while ($line=<in>) {
  chomp $line;

#	blastn -outfmt 6 ...
#	'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
#	comp9_c0_seq1	gi_109255272_ref_NC_008168.1__Choristoneura_occidentalis_granulovirus,_complete_genome	0.00	0	0	0	0	0	0	0	3e-19	99.0
  my @data = split /\t/, $line;

#  my ($front, $back) = split (/len[:]/, $data[0]);
#	there is no "len[:]" in $data[0] !!!
#  my @length = split (/[_]/, $back); 
#  my $length = $length[0];

#	who knows what the originally formated line was composed of
#	and that completely devalues all of the $data[...]

#	many of these columns are actually 0 so this will still fail!!!
  
#  if ($data[11]/$length >= 0.75) {
#  if ($data[11]/$data[3] >= 0.75) {
#	if data[3] is 0, skip?
  if (( $data[3] > 0 ) && ( $data[11]/$data[3] >= 0.75) ){
    $pline = "$data[0]\t$count{$data[0]}\t$data[1]\t$data[10]\t$data[11]\t$seq{$data[0]}";
    if (!$flag{$pline}) {
      print out "$pline\n";
      $flag{$pline} = 1;
    }
  }
}
close in;
close out;

