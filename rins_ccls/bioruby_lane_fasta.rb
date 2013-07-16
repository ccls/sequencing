#!/usr/bin/env ruby

require 'bio'	#	bioruby gem
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'sequencing_lib'

#	read a fasta file
#	find the reads that have both a left AND right lane (1 and 2)
#		(sometimes this is 0 and 1)
#	

if ARGV.empty?
	puts "#{$0} filename(s)\n\nAt least one input fasta file name is required.\n\n"
	exit
end

ARGV.each do |filename|
	inputfile     = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)
	lanes = [1,2]	#[]

	command = "grep '>' #{filename} | wc -l"
	puts "Counting sequences with ..."
	puts command
#	total_sequences = inputfile.count
	total_sequences = `#{command}`.to_i
#	will return string with leading spaces and trailing carriage return.  
#=> " 8899236\n"
#	to_i strips it all off and returns just the integer. Awesome.
##=> 8899236

	puts "Found #{total_sequences} sequences"

#	this is probably faster
#grep '>' ~/github_repo/ccls/sequencing/trinity_input_single.fasta  | wc -l
# 8899236
#	yep.  about 10-15x faster

	#	using +1 as log of 2-10 is 1, 11-100 is 2
	max_digits = Math.log10( total_sequences + 1 ).ceil

#	inputfile.each_with_index {|e,i|
#		printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", i+1, total_sequences
#		lane = e.definition.lane.to_i
#		lanes.push(lane) unless lanes.include?(lane)
#		sequences[e.definition.delane] ||= {}
#		#
#		#	This is saving the entire entry (name and sequence) which consumes 
#		#	more memory.  May have to only remember name as did before,
#		#	then rewind and reread file searching only for paired sequences.
#		#	
#		#		sequences[e.definition.delane] += 1
#		#
#		sequences[e.definition.delane][lane] = e
#
#
##
##Reading sequence 11319455 of 15421180
##
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta.rb:133: [BUG] object allocation during garbage collection phase
##ruby 2.0.0p195 (2013-05-14 revision 40734) [x86_64-linux]
##
##-- Control frame information -----------------------------------------------
##c:0011 p:---- s:0051 e:000050 CFUNC  :sub
##c:0010 p:0031 s:0046 e:000045 METHOD /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta.rb:133 [FINISH]
##c:0009 p:---- s:0042 e:000041 CFUNC  :new
##c:0008 p:0025 s:0038 e:000036 METHOD /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/splitter.rb:55
##c:0007 p:0085 s:0033 e:000032 METHOD /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile.rb:288
##c:0006 p:0020 s:0029 e:000028 METHOD /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile.rb:335 [FINISH]
##c:0005 p:---- s:0025 e:000024 CFUNC  :each_with_index
##c:0004 p:0130 s:0022 e:000021 BLOCK  /my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:42 [FINISH]
##c:0003 p:---- s:0007 e:000006 CFUNC  :each
##c:0002 p:0077 s:0004 E:000e08 EVAL   /my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:17 [FINISH]
##c:0001 p:0000 s:0002 E:0001e8 TOP    [FINISH]
##
##/my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:17:in `<main>'
##/my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:17:in `each'
##/my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:42:in `block in <main>'
##/my/home/jwendt/dna/bin/bioruby_lane_fasta.rb:42:in `each_with_index'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile.rb:335:in `each_entry'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile.rb:288:in `next_entry'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/splitter.rb:55:in `get_parsed_entry'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/splitter.rb:55:in `new'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta.rb:133:in `initialize'
##/my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta.rb:133:in `sub'
##
##-- C level backtrace information -------------------------------------------
##
##-- Other runtime information -----------------------------------------------
##
##* Loaded script: /my/home/jwendt/dna/bin/bioruby_lane_fasta.rb
##
##* Loaded features:
##
##    0 enumerator.so
##    1 /usr/lib64/ruby/enc/encdb.so
##    2 /usr/lib64/ruby/enc/trans/transdb.so
##    3 /usr/lib64/ruby/rbconfig.rb
##    4 /usr/share/rubygems/rubygems/compatibility.rb
##    5 /usr/share/rubygems/rubygems/defaults.rb
##    6 /usr/share/rubygems/rubygems/deprecate.rb
##    7 /usr/share/rubygems/rubygems/errors.rb
##    8 /usr/share/rubygems/rubygems/version.rb
##    9 /usr/share/rubygems/rubygems/requirement.rb
##   10 /usr/share/rubygems/rubygems/platform.rb
##   11 /usr/share/rubygems/rubygems/specification.rb
##   12 /usr/share/rubygems/rubygems/exceptions.rb
##   13 /usr/share/rubygems/rubygems/defaults/operating_system.rb
##   14 /usr/share/rubygems/rubygems/core_ext/kernel_gem.rb
##   15 /usr/share/rubygems/rubygems/core_ext/kernel_require.rb
##   16 /usr/share/rubygems/rubygems.rb
##   17 /usr/share/rubygems/rubygems/path_support.rb
##   18 /usr/share/rubygems/rubygems/dependency.rb
##   19 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio.rb
##   20 /my/home/jwendt/dna/bin/sequencing_lib.rb
##   21 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile.rb
##   22 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/buffer.rb
##   23 /usr/share/ruby/tsort.rb
##   24 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/autodetection.rb
##   25 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/common.rb
##   26 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/na.rb
##   27 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/aa.rb
##   28 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/compat.rb
##   29 /usr/share/ruby/cgi/util.rb
##   30 /usr/lib64/ruby/strscan.so
##   31 /usr/share/ruby/erb.rb
##   32 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/format.rb
##   33 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/sequence_masker.rb
##   34 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence.rb
##   35 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/reference.rb
##   36 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/feature.rb
##   37 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db.rb
##   38 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/sequence/dblink.rb
##   39 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta/defline.rb
##   40 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/fasta.rb
##   41 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/io/flatfile/splitter.rb
##   42 /usr/lib64/ruby/date_core.so
##   43 /usr/share/ruby/date/format.rb
##   44 /usr/share/ruby/date.rb
##   45 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/genbank/common.rb
##   46 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/genbank/genbank.rb
##   47 /my/home/jwendt/.gem/gems/bio-1.4.3.0001/lib/bio/db/genbank/genpept.rb
##
##[NOTE]
##You may have encountered a bug in the Ruby interpreter or extension libraries.
##Bug reports are welcome.
##For details: http://www.ruby-lang.org/bugreport.html
#
#
#	}
#	puts	#	mostly for newline after status

	#	This is blunt and effective, I think, but it all happens behind the scenes
	#	It would be nice to make this more verbose.
	#	It is substantially easier on memory though as it doesn't load
	#	the entire file into memory.
	puts "Scanning #{filename} for paired sequences (this can take a while) with ..."
	command = "grep '>' #{filename} | awk -F/ '{print $1}' | sort | uniq -d | sed 's/^>//'"
	puts command
	paired_sequences_a = `#{command}`.chomp.split
	paired_sequences = {}
	paired_sequences_a.each{|i| paired_sequences[i] = 1 }
	paired_sequences_a = nil



#	unpaired_sequences.sequence.select{|name,count| count == 1}
#	paired_sequences.sequence.select{  |name,count| count == 2}
#	confused_sequences.sequence.select{|name,count| count > 2 || count < 1 }
#
#	after this, most of this would have to change
#
#
#	unpaired_sequences = sequences.select{|k,v| v.length == 1 }
#	paired_sequences   = sequences.select{|k,v| v.length == 2 }
#	confused_sequences = sequences.select{|k,v| v.length > 2 || v.length < 1 }
#
#	puts "Found #{unpaired_sequences.length} unpaired sequences"

	puts "Found #{paired_sequences.length} paired sequences"

#	#	Not using them, so drop them to possibly free up so memory
#	unpaired_sequences = nil
#	confused_sequences = nil

	puts "Opening output files."
	streams = {}
	lanes.each do |lane|
		streams[lane] = File.open("#{root_filename}_#{lane}#{root_extname}",'w')
	end

	puts "Writing to output files."
	total_paired_sequences = paired_sequences.length
#	max_digits = Math.log10( total_paired_sequences + 1 ).ceil


#	May want to do normal for loop and pop values off to free up memory
#	if this is where there is a memory problem.


#	paired_sequences.each_with_index do |kv,i|
#		printf "\rWriting sequence %#{max_digits}d of %#{max_digits}d", i+1, total_paired_sequences
#		kv[1].each do |lane,sequence|
#			streams[lane].puts sequence
#		end
#	end
#	puts	#	for newline after status line

	inputfile.each_with_index do |entry,index|
#
#	This comparison is probably going to be SLOW for big files with many pairs.
#	Found 6693002 paired sequences
#	That's a big array to check over and over and over and over and ....
#
#
		printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", index+1, total_sequences
#	WAY WAY WAY TOO SLOW.  OMG
#		if paired_sequences.include?(entry.definition.delane)
#	Hash#has_key? IS WAY SUPER FASTER COMPARED TO Array#includes?
#	I would've thunk it the opposite.
		if paired_sequences.has_key?(entry.definition.delane)
#
#	Add some verbosity here?
#

			#	if matches, write to the correct file
			streams[entry.definition.lane.to_i].puts entry
#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#			streams[entry.definition.lane.to_i].puts entry
		end
	end
	puts
	puts "Closing output files."
	streams.each{|k,v|v.close}
end

__END__


This version is RIDICULOUSLY slow for big files.

ARGV.each do |filename|
	inputfile     = Bio::FlatFile.auto(filename)
	root_extname  = File.extname(filename)
	root_filename = File.basename(filename,root_extname)

	sequence_names = Hash.new(0)
	inputfile.each do |entry|
		#	HS1:209:c0v70acxx:6:1101:1633:2228_1:N:0:GATCAG
		#	=> HS1:209:c0v70acxx:6:1101:1633:2228
		#	HS1:209:c0v70acxx:6:1101:1633:2228 1:N:0:GATCAG
		#	=> HS1:209:c0v70acxx:6:1101:1633:2228
		#	HS2:360:D1NTUACXX:8:1101:12120:2042/1
		#	=> HS2:360:D1NTUACXX:8:1101:12120:2042
		#	comp2_c0_seq1 len=295 path=[296:0-50 101:51-172 347:173-294]
		#	=> comp2_c0_seq1 len=295 path=[296:0-50 101:51-172 347:173-294]
		#	
		sequence_names[entry.definition.delane] += 1
	end

	#	only those that were there twice (left and right)
	sequence_names.delete_if{|name,count| count < 2}

	#	convert to array of just the names
	sequence_names = sequence_names.collect{|name,count| name }

	#	open the output files and keep in order in this array
	streams = [File.open("#{root_filename}_1#{root_extname}",'w'),
		File.open("#{root_filename}_2#{root_extname}",'w')]

#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#	streams = [File.open("#{root_filename}_0#{root_extname}",'w'),
#		File.open("#{root_filename}_1#{root_extname}",'w'),
#		File.open("#{root_filename}_2#{root_extname}",'w')]

	inputfile.rewind	#	already read the file, so need to rewind
	inputfile.each do |entry|
		if sequence_names.include?(entry.definition.delane)
			#	if matches, write to the correct file
			streams[entry.definition.lane.to_i - 1].puts entry
#	if 0/1 become actual possiblities, open 3 streams?
#	should end up with 1 empty file at the end
#			streams[entry.definition.lane.to_i].puts entry
		end
	end
	streams.each{|f|f.close}

end	#	ARGV.each do |filename|
__END__
