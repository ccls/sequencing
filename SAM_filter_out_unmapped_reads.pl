#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;



main: {

	my $sam_reader = new SAM_reader($sam_file);
    
    
    my $filtered_count = 0;
    my $total_count = 0;


#print STDERR "JAKE - BEGIN\n";
#print STDERR "JAKE - sam file: " . $sam_file . "\n";
#print STDERR "JAKE - file size: " . ((-s $sam_file)||"") . "\n";

#	my printed output needs to go to STDERR
#	as the calling script redirects the output of this to a file. 

	while ($sam_reader->has_next()) {

#print STDERR "JAKE - in has next while loop.\n";
		
		my $sam_entry = $sam_reader->get_next();
        $total_count++;
        
        if ($sam_entry->is_query_unmapped()) {
            $filtered_count++;
        }
        else {
#print STDERR "SAM ENTRY: " . $sam_entry->toString() . "\n";
            print $sam_entry->toString() . "\n";
        }
    }

#	When running the 3 Fallon datasets and either of the CAHPV, all will fail with
#	something like ... on the last line of the file.
#	Illegal division by zero at /Users/jakewendt/rins/trinity/util/../util/SAM_filter_out_unmapped_reads.pl line 37, <$__ANONIO__> line 300639.
#	
#	The calling script actually has the path to this script hard coded
#	so this must be modified in place or linked.
#	/Users/jakewendt/rins/trinity/util/SAM_filter_out_unmapped_reads.pl

#-filtered 155 of 324 reads as unaligned = 47.84% unaligned reads

#print STDERR "JAKE - total count:" . $total_count . "\n";
if( $total_count > 0 ){
    print STDERR "-filtered $filtered_count of $total_count reads as unaligned = " . sprintf("%.2f", $filtered_count / $total_count * 100) . "\% unaligned reads\n";
}
    
#print STDERR "JAKE - END\n";

	exit(0);

}
