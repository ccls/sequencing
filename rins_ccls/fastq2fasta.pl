#!/usr/bin/env perl
##!/usr/bin/perl

# This script is used to change .fastq file to .fasta file

#
#	CCLS Modified
#	
#	This script REQUIRES that the sequence is on a single line.
#	Not sure if that is standard, but always seems to be the case.
#	
#	The original version did not remove the @ symbol from the
#	sequence name which causes all sorts of problems if the data
#	makes it into a SAM format.
#
#	Most of these scripts also assume the Illumina sequence name
#	format with a trailing /1 or /2 for lane identification.
#
#	Our data appears to be in Casava format and this script
#	replaces the SPACE field separator with an underscore.
#	Not sure if this is an issue, but some other scripts do
#	seem to be written to expect the /1 or /2.
#
#	Perhaps fix that?
#

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

if ($#ARGV ne 1) {
	print "command line example: perl fastq2fasta.pl input.fastq output.fasta\n";
	exit;
}

$infile = $ARGV[0];
$outfile = $ARGV[1];


open (in, "<$infile");
open (out, ">$outfile");
while ($line=<in>) {
 
	$count ++;
	if ($count%1000000 eq 0) {
		print "$count\n";
	}

	chomp $line;
	$line=~s/\s+/\_/g;  #	sequence name (replace spaces with underscores. Why?)
	$line=~s/^@//;      #	sequence name (remove that leading @ symbol)
	print out ">$line\n";

	$line = <in>;	#	sequence
	print out $line;

	$line = <in>;	#	ignore 3rd FASTQ line
	$line = <in>;	#	ignore 4th FASTQ line
}

close in;
close out;
