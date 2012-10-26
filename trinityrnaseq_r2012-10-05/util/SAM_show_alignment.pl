#!/usr/bin/env perl 

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use SAM_reader;
use SAM_entry;


my $FASTA_LENGTH = 60;

my $usage = "usage: $0 alignments.sam target.fasta [JUST_ALIGN_STATS=0]\n\n";

my $sam_alignments = $ARGV[0] or die $usage;
my $target_fasta = $ARGV[1] or die $usage;
my $JUST_STATS = $ARGV[2] || 0;

main: {

    my $fasta_reader = new Fasta_reader($target_fasta);
    my %fasta_seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    my $sam_reader = new SAM_reader($sam_alignments);

    my $curr_scaff_acc = "";
    my $curr_scaff_seq = "";

    while (my $sam_entry = $sam_reader->get_next() ) {
        
        my $read_name = $sam_entry->reconstruct_full_read_name();
        my $read_seq = uc $sam_entry->get_sequence();
        my $scaff_name = $sam_entry->get_scaffold_name();
        
        if ($scaff_name ne $curr_scaff_acc) {
            $curr_scaff_acc = $scaff_name;
            $curr_scaff_seq = uc $fasta_seqs{$scaff_name};
        }
        
        my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
        
        unless (scalar(@$genome_coords_aref) >= 1) { next; }
        
        print "// $scaff_name (top) vs. $read_name (bottom):\n\n" unless $JUST_STATS;
        &draw_alignment($scaff_name, $read_name, $genome_coords_aref, $read_coords_aref, \$curr_scaff_seq, \$read_seq);
        
        
    }


    exit(0);
}



####
sub draw_alignment {
    my ($scaff_name, $read_name, $genome_coords_aref, $read_coords_aref, $curr_scaff_seq_sref, $read_seq_sref) = @_;
    
    my $scaff_align_string = "";
    my $read_align_string = "";
    
    my $prev_scaff_coord = 0;
    my $prev_read_coord = 0;

    my $num_indels = scalar(@$genome_coords_aref) - 1;
    my $num_aligned_bases = 0;

    while (@$genome_coords_aref) {
        my $genome_coordset = shift @$genome_coords_aref;
        my $read_coordset = shift @$read_coords_aref;

        my ($genome_lend, $genome_rend) = @$genome_coordset;
        my ($read_lend, $read_rend) = @$read_coordset;

        $num_aligned_bases += $read_rend - $read_lend + 1;

        if ($prev_read_coord > 0 && $read_lend - $prev_read_coord > 1) {
            ## add gaps to read string
            for (my $i = $prev_read_coord + 1; $i < $read_lend; $i++) {
                $read_align_string .= substr($$read_seq_sref, $i-1, 1);
                $scaff_align_string .= "-";;
            }
        }
        
        if ($prev_scaff_coord > 0 && $genome_lend - $prev_scaff_coord > 1) {
            # add gaps to genome string
            for (my $i = $prev_scaff_coord + 1; $i < $genome_lend; $i++) {
                $scaff_align_string .= substr($$curr_scaff_seq_sref, $i-1, 1);
                $read_align_string .= "-";
            }
        }
        
        $prev_read_coord = $read_rend;
        $prev_scaff_coord = $genome_rend;
        
        my $read_seq_region = substr($$read_seq_sref, $read_lend - 1, $read_rend - $read_lend + 1);
        my $genome_seq_region = substr($$curr_scaff_seq_sref, $genome_lend - 1, $genome_rend - $genome_lend + 1);
        
        $scaff_align_string .= $genome_seq_region;
        $read_align_string .= $read_seq_region;
    
        
    
    }

    my $alignment_length = length($scaff_align_string);
    my $num_gaps = 0;
    while ($scaff_align_string =~ /-/g) { $num_gaps++; }
    while ($read_align_string =~ /-/g) { $num_gaps++; }

    my $percent_gap = $num_gaps / $alignment_length * 100;
        

    my ($num_mismatches) = &print_pretty($scaff_align_string, $read_align_string);
    
    print join("\t", "#", "scaff_name", "read_name", "aligned_bases", "matches", "mismatches", "indel_bkpts", "sum_indel_lens", "pct_mismatches", "pct_indel_bkpts", "pct_indel_lens") . "\n";
    print join("\t", "#", $scaff_name, $read_name,  
               $num_aligned_bases, $num_aligned_bases - $num_mismatches, $num_mismatches, $num_indels, $num_gaps, 
               sprintf("%.2f", $num_mismatches/$num_aligned_bases*100),
               sprintf("%.2f", $num_indels/$num_aligned_bases*100),
               sprintf("%.2f", $percent_gap),
               ) . "\n\n\n";
    
    return;
    
    
}
        
####
sub print_pretty {
    my ($scaff_align_string, $read_align_string) = @_;

    if (length($scaff_align_string) != length($read_align_string) ) {
        die "Error, alignment lengths differ.";
    }
    
    my @scaff_align_chars = split(//, $scaff_align_string);
    my @read_align_chars = split(//, $read_align_string);
    
    
    my $top_text = "";
    my $match_text = "";
    my $bottom_text = "";

    my $num_mismatches = 0;

    for (my $i = 0; $i <= $#scaff_align_chars; $i++) {
        
        if ($i != 0 && $i % 60 == 0) {
            
            print join("\n", $top_text, $match_text, $bottom_text) . "\n\n" unless $JUST_STATS;
            
            $top_text = ""; $match_text = ""; $bottom_text = "";
    
        }

        my $match_char = ".";
        if ($scaff_align_chars[$i] ne "-" && $read_align_chars[$i] ne "-") {

            $match_char = ($scaff_align_chars[$i] eq $read_align_chars[$i]) ? "|" : "*";

            if ($match_char eq "*") {
                $num_mismatches++;
            }
            
        }
        $top_text .= $scaff_align_chars[$i];
        $bottom_text .= $read_align_chars[$i];
        $match_text .= $match_char;

    }
    
    if ($top_text) {
        print join("\n", $top_text, $match_text, $bottom_text) . "\n\n" unless $JUST_STATS;
    }


    return ($num_mismatches);
    
}
