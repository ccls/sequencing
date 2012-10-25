#!/usr/bin/env perl
##!/usr/bin/perl
#
# This script is used to compress reads used the LZW compression ratio as a hard cutoff
#
# Use: input file, direct to output with > 

use warnings;
use strict;

my $input = $ARGV[0];
my $bits = 8;
my $compress_ratio_thrd = $ARGV[1];	#	This is usually 0.5

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
#
#	The above is removing the trailing "lane specifier", but it is
#	not standard.  Our lane data is different.
#
#	(from my ruby script)
#		# >@HWI-ST281_0133:3:1:1222:2139#0/1
#		name.gsub!(/^>/,'')
#		# @HWI-ST281_0133:3:1:1222:2139#0/1
#		# 
#		# @HWI-ST281_0133:3:1:1222:2139#0/1
#		name.gsub!(/\/\d+$/,'')
#		# @HWI-ST281_0133:3:1:1222:2139#0
#		# 
#		# @HWI-ST977:132:C09W8ACXX:7:2307:10304:95858_2:N:0:CGTAGG
#		name.gsub!(/_\d{1}:.*$/,'')
#		# @HWI-ST977:132:C09W8ACXX:7:2307:10304:95858
#		# 
#
#	Adding this to "delane" our data...
		$namesfastq =~s/_\d{1}:.*$//;

#
#	I have modified the scripts that read this output
#	to "delane" it before doing anything
#	so this may not be entirely necessary now.
#


		$namesfastq =~s/^\>//;

		print "$namesfastq\n";
	}
	@dictionary = ();
	%dictionary = ();

}


close (IN);
