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
	:output_filename          => 'results.txt',
	:output_suffix            => 'dark',
	:bowtie_version           => 1,
	:link_sample_fa_files     => false,
	:pre_chopped              => false,
	:chop_read_length         => 25,
	:minIdentity              => 80,
	:compress_ratio_thrd      => 0.5,
	:iteration                => 2,
	:bowtie_threads           => 6,
	:bowtie_mismatch          => 3,
	:paired_fragment_length   => 300,
	:min_contig_length        => 300,
	:trinity_threads          => 6,
	:blastn_evalue_thrd       => 0.05,
	:similarity_thrd          => 0.8,
	:mailto                   => '',
	:die_on_failed_file_check => false,
	:files                    => {}
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
EOB
		exit
	end
end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options. What's left is the list of files to resize.
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
		outbase = ''	#	outside to save last value for pull_reads_fasta

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

			command = "bowtie2 -p 4 -N 1 " <<
				"#{(file_format == 'fastq')? '-q' : '-f'} " <<
				"-x #{bowtie_human_index} " <<
				"-1 #{prevbase}.1.#{file_format} " <<
				"-2 #{prevbase}.2.#{file_format} " <<
				"-S #{outbase}.sam " <<
				"--all " <<
				"--un-conc #{outbase}.#{file_format}"
			command.execute

			"#{outbase}.sam".file_check(die_on_failed_file_check)
			"#{outbase}.1.#{file_format}".file_check(die_on_failed_file_check)
			"#{outbase}.2.#{file_format}".file_check(die_on_failed_file_check)
	
		end	#	%w( hg18 hg19 ).each do |hg|


#	mapper_repeatmasker.py
#		mega1="/root/blast-2.2.23/bin/megablast -i /mnt/mfiles/reads/" + fFASTA	
#		mega1=mega1 + " -d \"/mnt/mfiles/ref/"
#		mega1=mega1 + database_name
#		mega1=mega1 + "\" -a 4 -D 2 -e 0.0000001 -W 16 -b 0 -v 5 -f -F F -o /mnt/mfiles/reads/"
#		mega1=mega1 + fOUT


#	mapper_postunmapped.py
#	#Megablast on reads over Bacterial Genome#############RAM Tools####################### 
#	mega1="/root/blast-2.2.23/bin/megablast -i /mnt/mfiles/reads/" + fFASTA	
#	mega1=mega1 + " -d \"/mnt/mfiles/ref/bacterial.db\" -m 7 -a 4 -D 2 -e 0.0000001 -W 16 -b 0 -v 5 -f -F F -o /mnt/mfiles/reads/"
#	mega1=mega1 + fOUT

#	mapper_postvelvet.py	
#	mega1="/root/blast-2.2.23/bin/megablast -i /mnt/mfiles/reads/" + fFASTA	
#	mega1=mega1 + " -d \"/mnt/mfiles/ref/bacterial.db\" -m 7 -a 4 -D 2 -e 0.0000001 -W 16 -b 0 -v 5 -f -F F -o /mnt/mfiles/reads/"
#	mega1=mega1 + fOUT

#	all ...
#		-a 4
#			-num_threads 4	#	This seems to always raise 
#    "/Users/jakewendt/sequencing/ncbi-blast-2.2.27+-src/c++/src/corelib/ncbiobj.cpp", 
#			line 689: 
#			Critical: ncbi::CObject::ThrowNullPointerException() - Attempt to access NULL pointer.
#		-D 2
#			-outfmt ???? kinda
#		-e 0.0000001
#			-evalue 0.0000001
#		-W 16
#			-word_size 16
#		-b 0
#			-num_alignments 0		( not with outfmt 6 ??? )
#		-v 5
#			-num_descriptions 5
#		-f ( show full ids in the output )
#			
#		-F F ( Don't filter query sequence )
#			-dust no  ????			

#		-m 7
#			-outfmt 5	( xml output )	#	causes ...
#			Warning: The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output
#			BLAST query/options error: No hits are being saved
#			Using max_target_seqs causes ...
#			Error: Argument "num_alignments". Incompatible with argument:  `max_target_seqs'
#		Which defeats the purpose? num_alignments=0 doesn't just report the
#		unaligned so I really don't get this at all. So confused.


#	Really confused. This will produce a file that will require heavy parsing. 
#	In addition, doesn't just include the misses.  Lots of hits to database that
#	has already said that it didn't hit?  I'm guess that it is because it is not
#	pair checking? Need to read the Parse???.cc file and see what PathSeq keeps.

#	s = The blastn output (freaking huge)
#	> blastn -db /Volumes/cube/working/indexes/Homo_sapiens.GRCh37.69.cdna.all -query trinity_output/both.fa -evalue 0.0000001 -word_size 16 -num_alignments 0 -outfmt 0 -num_descriptions 5 > mega_blast_both_fa.out &
#
#	s.scan(/Query= (.*)\n\nLength=\d+\n\n\n.*No hits/)
#	=> [["HWI-ST977:132:C09W8ACXX:7:1102:6838:8469/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:17366:8499/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:19346:8490/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:2558:8566/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:14765:8297/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:1861:8665/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:18376:8345/1"], ["HWI-ST977:132:C09W8ACXX:7:1102:4826:8709/1"]]

#
#	The problem with this is that we would read the entire file.
#	Probably better to read line by line, save the sequence name from 
#		"Query= HWI-ST977:132:C09W8ACXX:7:1102:6838:8469/1"
#	And then use it if find 
#		"*****  NO hits found *****"
#	After that, delane and pull_reads_from() ...
#

#	OK, this blastn has taken at least 48 hours.  Too long to be useful.

#	used legacy_blast.pl to convert to new command (using multiple threads causes crash)
#	legacy_blast.pl megablast -i trinity_output/both.fa -d /Volumes/cube/working/indexes/Homo_sapiens.GRCh37.69.cdna.all -a 4 -D 2 -e 0.0000001 -W 16 -b 0 -v 5 -f -F F -o somefile --path /Users/jakewendt/RINS_BASE/bin --print_only
#	blastn -db /Volumes/cube/working/indexes/Homo_sapiens.GRCh37.69.cdna.all -query trinity_output/both.fa -evalue 0.0000001 -word_size 16 -num_alignments 0 -outfmt 5 -max_target_seqs 5 > mega_blast_both_fa.xml


		
#		pull_reads_from_fastas(
#			files.keys.sort.collect{|k| "#{k}#{outbase}.names" },
#			files.keys.sort.collect{|k| "#{k}lane.fa" },
#			files.keys.sort.collect{|k| "#{k}lane_bowtie_non_human.fasta" })
#
#	link last files as 
#	files.keys.sort.collect{|k| "#{k}lane_bowtie_non_human.fasta" })

#		files.keys.sort.each_with_index{|k,i| 
#			FileUtils.ln_s("#{outbase}.#{i+1}.#{file_format}","#{k}lane_bowtie_non_human.#{file_format}") }

		puts "de novo assembly using Trinity"
#		command = "Trinity.pl --seqType fa " <<
		command = "Trinity.pl --seqType #{(file_format == 'fastq')? 'fq' : 'fa'} " <<
			"--group_pairs_distance #{paired_fragment_length} " <<
			"--min_contig_length #{min_contig_length} " <<
			"--output trinity_output " <<
			"--CPU #{trinity_threads} " <<
			"--bfly_opts \"--stderr\" --JM 1G "
		files.keys.sort.each_with_index{|k,i| 
			command << "--#{k} #{outbase}.#{i+1}.#{file_format} " }

#		files.each_pair { |k,v| command << "--#{k} #{k}lane_bowtie_non_human.#{file_format} " }


#
#	will fastq in produce fastq out with trinity????
#

		command.execute
		"trinity_output/Trinity.fasta".file_check(die_on_failed_file_check)
		FileUtils.cp("trinity_output/Trinity.fasta","trinity_non_human.fasta")

	end

end


#########

system "date";
darkness = Darkness.new(o)
darkness.prepare_output_dir_and_log_file


#darkness.file_format_check_and_conversion


darkness.bowtie_non_human
darkness.blastn_non_human("trinity_non_human.fasta")
darkness.wrap_things_up





__END__
