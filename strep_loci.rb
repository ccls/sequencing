#!/usr/bin/env ruby

require 'erb'
require 'yaml'
require 'optparse'
require 'fileutils'

o = {
	:in_dir  => "/Volumes/cube/working/samples/targeted",
	:out_dir => "./"
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	#	on -h -help --help
	#	This will be followed by the command line options.
	opts.banner =<<EOB
Usage: #{File.basename($0)} [options]
EOB

	# Define the options, and what they do
	opts.on( '-i','--in DIR_NAME', 'input dir') do |s|
		o[:in_dir] = s
	end

	opts.on( '-o','--out DIR_NAME', 'output dir') do |s|
		o[:out_dir] = s
	end

	opts.on( '-f','--file FILE_BASE_NAME', 'in files') do |s|
		o[:in_file] = s
	end

	opts.on( '-h', 'Display simple help screen' ) do
		puts opts
		puts
		exit
	end
end
 
# Parse the command-line. Remember there are two forms
# of the parse method. The 'parse' method simply parses
# ARGV, while the 'parse!' method parses ARGV and removes
# any options found there, as well as any parameters for
# the options. What's left is the list of files to resize.
optparse.parse!

#	file required
unless o[:in_file]
	puts optparse	#	Basically display the command line help
	exit
end




#1run create alignment using bowtie to strep loci reference
command = "bowtie2 -a -q -N 0 --qc-filter " <<
	"-x /Volumes/cube/working/indexes/strep_loci " <<
	"-1 #{o[:in_dir]}/#{o[:in_file]}_R1.fastq " <<
	"-2 #{o[:in_dir]}/#{o[:in_file]}_R2.fastq " <<
	"-S #{o[:out_dir]}/#{o[:in_file]}_N0.sam"
system(command)


#2Convert Sam to Bam file
command = "samtools view -uS #{o[:out_dir]}/#{o[:in_file]}_N0.sam " << 
	"| samtools sort - #{o[:out_dir]}/#{o[:in_file]}_N0"
system(command)

#3 create index from bam file
command = "samtools index #{o[:out_dir]}/#{o[:in_file]}_N0.bam"
system(command)



__END__

It would be great if I could give it a folder of input files (fastq with same name with  _R1 and _R2 designation) and an output folder to dump the output into.

Also, If I could change the Bowtie2 settings (options and index) easily that would be great. Thank you so much for the help, this will take me hours/days.... I doubt it will take you a half hour.


Btw... I'm finally feeling better. I got the office cold. 
