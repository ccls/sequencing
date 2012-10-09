#!/usr/bin/env ruby

class FileFormatDetector
	attr_accessor :filename, :format
	def initialize(filename)
		self.filename = filename
		self.format   = detect_format
	end
	def detect_format
		self.format = if File.exists?(filename)
			if File.file?(filename)
				File.open(filename, 'r') do |f|
					line = f.gets
					if( line =~ /^@/ )
						'fastq'
					elsif( line =~ /^>/ )
						'fasta'
					else
						'unknown'
					end
				end
			else
				'notafile'
			end
		else
			'nonexistant'
		end
	end
end


ARGV.each do |infilename|
	file_format = FileFormatDetector.new(infilename)
	puts "#{file_format.filename}:#{file_format.format}"
end
