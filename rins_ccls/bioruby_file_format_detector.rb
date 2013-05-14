#!/usr/bin/env ruby

require 'bio'	#	bioruby gem

if ARGV.empty?
	puts "#{$0} filename(s)"
	puts
	puts "At least one file name is required."
	puts
	exit
end

ARGV.each do |filename|

	inputfile = Bio::FlatFile.auto(filename)

	format = case inputfile.dbclass.to_s
		when 'Bio::FastaFormat' then 'fasta'
		when 'Bio::Fastq' then 'fastq'
		else 'unknown'
	end
	
	puts "#{filename}:#{format}"

end	#	ARGV.each do |filename|

__END__


This is so simple, it may not even be worthy of a script

