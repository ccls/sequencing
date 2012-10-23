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
					while( line = f.gets )
						#	skip over comment lines
						next if( line =~ /^#/ ) #	or whatever comment character is.

						type = if( line =~ /^@/ )
							'fastq'
						elsif( line =~ /^>/ )
							'fasta'
						else
							'unknown'
						end

						#	If made it this far, we're done.
						#	Apparently break will return a value.
						break	type 
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

#
#	This file can be loaded as a library, but when used
#	at the command line it will check each ARGV.
#
if __FILE__==$0

	usage =<<EOUSAGE

Usage: #{File.basename($0)} filenames

Checks file's first character to determine if fasta or fastq.

EOUSAGE

	#	show help if arg is empty or 
	#	used -h or --help or even --hotdog
	if( ( ARGV.empty? ) or
		( !ARGV.empty? and ARGV[0].match(/-h/i) ) )
		puts usage
		exit
	end

	ARGV.each do |infilename|
		file_format = FileFormatDetector.new(infilename)
		puts "#{file_format.filename}:#{file_format.format}"
	end
end
