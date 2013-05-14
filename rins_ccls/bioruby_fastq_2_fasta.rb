#!/usr/bin/env ruby

require 'bio'

if ARGV.empty?
	puts "#{$0} input_fastq <output_fasta>\n\nAt least input file name is required.\n\n"
	exit
end

inputfile = Bio::FlatFile.auto( ARGV.shift )

if inputfile.dbclass.to_s != 'Bio::Fastq'
	puts "\nInput file is not a fastq\n\n"
	exit
end

o = ( ARGV.empty? ) ? STDOUT : File.open(ARGV.shift,'w')

inputfile.each { |e| o.puts e.seq.to_fasta(e.definition) }

o.close unless o == STDOUT
