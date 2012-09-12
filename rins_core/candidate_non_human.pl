#!/usr/bin/perl
#
# used to find candidate non human species and their counts
#
###########################################################

my %count = ();


for ($i=0; $i<=$#ARGV; $i++) {
  open (in, "<$ARGV[$i]");
  while ($line=<in>) {
    chomp $line;
    $name = $line;
    $flag{$name} = 1;
  }
  close in;
}



if ($#ARGV eq 1) {

  %counted = ();
  open (in, "<chopped_leftlane.psl");
  for ($i=1; $i<=4; $i++) {
    $line = <in>;
  }
  while ($line=<in>) {
    chomp $line;
    @data = split /\t/, $line;
    $name = $data[9];
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;
    
    $species_name = $data[13];
    
    if ($flag{$name}) {
      if (!$counted{$species_name}{$name}) {
	$counted{$species_name}{$name} = 1;
	$count{$species_name} ++;
      }
    }
  }
  close in;
  
  %counted = ();
  open (in, "<chopped_rightlane.psl");
  for ($i=1; $i<=4; $i++) {
    $line = <in>;
  }
  while ($line=<in>) {
    chomp $line;
    @data = split /\t/, $line;
    $name = $data[9];
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;
    
    $species_name = $data[13];
    
    if ($flag{$name}) {
      if (!$counted{$species_name}{$name}) {
	$counted{$species_name}{$name} = 1;
	$count{$species_name} ++;
      }
    }
  }
  close in;
}

else {
  
  %counted = ();
  open (in, "<chopped_singlelane.psl");
  for ($i=1; $i<=4; $i++) {
    $line = <in>;
  }
  while ($line=<in>) {
    chomp $line;
    @data = split /\t/, $line;
    $name = $data[9];
    $name =~s/\/[0-9]+$//;
    $name =~s/^\>//;
    
    $species_name = $data[13];
    
    if ($flag{$name}) {
      if (!$counted{$species_name}{$name}) {
	$counted{$species_name}{$name} = 1;
	$count{$species_name} ++;
      }
    }
  }
  close in;
  
}

open (out, ">candidate_non_human.txt");
print out "non_human_species\traw_read_counts\n";
@species_names = sort {$count{$b} <=> $count{$a}} keys %count;
foreach $species_name (@species_names) {
  print out "$species_name\t$count{$species_name}\n";
}
close out;
  
