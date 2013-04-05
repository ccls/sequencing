#!/usr/bin/env ruby

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'
require 'ostruct'

$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'


puts "Running #{$0}"
o = {
	:config_filename          => 'config.yml',
#	:output_filename          => 'results.txt',
	:output_suffix            => 'dark',
#	:bowtie_version           => 1,
#	:link_sample_fa_files     => false,
#	:pre_chopped              => false,
#	:chop_read_length         => 25,
#	:minIdentity              => 80,
#	:compress_ratio_thrd      => 0.5,
#	:iteration                => 2,
#	:bowtie_threads           => 6,
#	:bowtie_mismatch          => 3,
#	:paired_fragment_length   => 300,
#	:min_contig_length        => 300,
#	:trinity_threads          => 6,
#	:blastn_evalue_thrd       => 0.05,
#	:similarity_thrd          => 0.8,
#	:mailto                   => '',
#	:die_on_failed_file_check => false,
#	:files                    => {}
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB
Usage: #{File.basename($0)} [options]
EOB

	# Define the options, and what they do

	opts.on( '--output_suffix STRING', 
		"Output directory suffix (default '#{o[:output_suffix]}')" ) do |s|
		o[:output_suffix] = s
	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end

	opts.on( '--help', 'Display thorough help screen' ) do
		puts opts
puts <<EOB

Usage: #{File.basename($0)} --output_suffix "darkdev.fastq.bowtie2.trinity20120608"

Loops through all of the bowtie_human_indexes in the config.yml
stripping out all of the matches using bowtie2's --un-conc option.

Trinity attempts to assemble the resultant file.

Those results are then blastn'd.

EOB
		exit
	end
end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options.
optparse.parse!

config = YAML::load( ERB.new( IO.read( File.join( o[:config_filename] ) ) ).result)
o.update( config )


#
#
#		Now we begin ....
#
#

class Darkness < CclsSequencer

#
#	This was designed to work with singles and pairs.
#	Now it really kinda requires pairs.
#

	def bowtie_non_human
		outbase = ''	#	outside to save last value for pull_reads_fasta ( irrelevant now )

		outbase = 'raw'
		files.keys.sort.each_with_index{|k,i| 
			FileUtils.ln_s(files[k],"raw.#{i+1}.#{file_format}") }
		"raw.1.#{file_format}".file_check(die_on_failed_file_check)
		"raw.2.#{file_format}".file_check(die_on_failed_file_check)

		[bowtie_human_indexes].flatten.each_with_index do |bowtie_human_index,i|
			prevbase = String.new(outbase)	#	DO NOT JUST SET = AS WILL NOT CREATE NEW STRING!
			name = File.basename(bowtie_human_index)
			outbase << "_not" if i == 0
			outbase << "_#{name}"
#	
#	Bowtie 2 version 2.0.0-beta7 by Ben Langmead (blangmea@jhsph.edu)
#	Usage: 
#	  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
#	
#	  <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
#	             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
#	  <m1>       Files with #1 mates, paired with files in <m2>.
#	             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#	  <m2>       Files with #2 mates, paired with files in <m1>.
#	             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#	  <r>        Files with unpaired reads.
#	             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#	  <sam>      File for SAM output (default: stdout)
#	
#	  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
#	  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
#	
#	Options (defaults in parentheses):
#	
#	 Input:
#	  -q                 query input files are FASTQ .fq/.fastq (default)
#	  --qseq             query input files are in Illumina's qseq format
#	  -f                 query input files are (multi-)FASTA .fa/.mfa
#	  -r                 query input files are raw one-sequence-per-line
#	  -c                 <m1>, <m2>, <r> are sequences themselves, not files
#	  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
#	  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
#	  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
#	  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
#	  --phred33          qualities are Phred+33 (default)
#	  --phred64          qualities are Phred+64
#	  --int-quals        qualities encoded as space-delimited integers
#	
#	 Presets:                 Same as:
#	  For --end-to-end:
#	   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
#	   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
#	   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
#	   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#	
#	  For --local:
#	   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#	   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#	   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
#	   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#	
#	 Alignment:
#	  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
#	  -L <int>           length of seed substrings; must be >3, <32 (22)
#	  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
#	  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
#	  --dpad <int>       include <int> extra ref chars on sides of DP table (15)
#	  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
#	  --ignore-quals     treat all quality values as 30 on Phred scale (off)
#	  --nofw             do not align forward (original) version of read (off)
#	  --norc             do not align reverse-complement version of read (off)
#	
#	  --end-to-end       entire read must align; no clipping (on)
#	   OR
#	  --local            local alignment; ends might be soft clipped (off)
#	
#	 Scoring:
#	  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
#	  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
#	  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
#	  --rdg <int>,<int>  read gap open, extend penalties (5,3)
#	  --rfg <int>,<int>  reference gap open, extend penalties (5,3)
#	  --score-min <func> min acceptable alignment score w/r/t read length
#	                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)
#	
#	 Reporting:
#	  (default)          look for multiple alignments, report best, with MAPQ
#	   OR
#	  -k <int>           report up to <int> alns per read; MAPQ not meaningful
#	   OR
#	  -a/--all           report all alignments; very slow, MAPQ not meaningful
#	
#	 Effort:
#	  -D <int>           give up extending after <int> failed extends in a row (15)
#	  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
#	
#	 Paired-end:
#	  -I/--minins <int>  minimum fragment length (0)
#	  -X/--maxins <int>  maximum fragment length (500)
#	  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
#	  --no-mixed         suppress unpaired alignments for paired reads
#	  --no-discordant    suppress discordant alignments for paired reads
#	  --no-dovetail      not concordant when mates extend past each other
#	  --no-contain       not concordant when one mate alignment contains other
#	  --no-overlap       not concordant when mates overlap at all
#	
#	 Output:
#	  -t/--time          print wall-clock time taken by search phases
#	  --un <path>           write unpaired reads that didn't align to <path>
#	  --al <path>           write unpaired reads that aligned at least once to <path>
#	  --un-conc <path>      write pairs that didn't align concordantly to <path>
#	  --al-conc <path>      write pairs that aligned concordantly at least once to <path>
#	  (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
#	  --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
#	  --quiet            print nothing to stderr except serious errors
#	  --met-file <path>  send metrics to file at <path> (off)
#	  --met-stderr       send metrics to stderr (off)
#	  --met <int>        report internal counters & metrics every <int> secs (1)
#	  --no-head          supppress header lines, i.e. lines starting with @
#	  --no-sq            supppress @SQ header lines
#	  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
#	  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
#	                     Note: @RG line only printed when --rg-id is set.
#	
#	 Performance:
#	  -o/--offrate <int> override offrate of index; must be >= index's offrate
#	  -p/--threads <int> number of alignment threads to launch (1)
#	  --reorder          force SAM output order to match order of input reads
#	  --mm               use memory-mapped I/O for index; many 'bowtie's can share
#	
#	 Other:
#	  --qc-filter        filter out reads that are bad according to QSEQ filter
#	  --seed <int>       seed for random number generator (0)
#	  --version          print version information and quit
#	  -h/--help          print this usage message

			command = "bowtie2 -N 1 " <<
				"#{(file_format == 'fastq')? '-q' : '-f'} " <<
				"-x #{bowtie_human_index} " <<
				"-1 #{prevbase}.1.#{file_format} " <<
				"-2 #{prevbase}.2.#{file_format} " <<
				"-S /dev/null " <<
				"--threads #{bowtie_threads} " <<
				"--un-conc #{outbase}.#{file_format}"
			command.execute
#				"-S #{outbase}.sam " <<

#			"#{outbase}.sam".file_check(die_on_failed_file_check)
			"#{outbase}.1.#{file_format}".file_check(die_on_failed_file_check)
			"#{outbase}.2.#{file_format}".file_check(die_on_failed_file_check)
	
		end	#	%w( hg18 hg19 ).each do |hg|

		puts "de novo assembly using Trinity"
#		command = "Trinity.pl --seqType fa " <<
#			"--group_pairs_distance #{paired_fragment_length} " <<
#			"--min_contig_length #{min_contig_length} " <<
#			"--CPU #{trinity_threads} " <<
#			"--bfly_opts \"--stderr\" --JM 1G "
		command = "Trinity.pl --seqType #{(file_format == 'fastq')? 'fq' : 'fa'} " <<
			"--output trinity_output " <<
			"--JM 2G "
		files.keys.sort.each_with_index{|k,i| 
			command << "--#{k} #{outbase}.#{i+1}.#{file_format} " }

		command.execute
		"trinity_output/Trinity.fasta".file_check(die_on_failed_file_check)
		FileUtils.cp("trinity_output/Trinity.fasta","trinity_non_human.fasta")

	end

end


#########

system "date";
darkness = Darkness.new(o)
darkness.prepare_output_dir_and_log_file
darkness.bowtie_non_human
darkness.blastn_non_human("trinity_non_human.fasta")
darkness.wrap_things_up

__END__
