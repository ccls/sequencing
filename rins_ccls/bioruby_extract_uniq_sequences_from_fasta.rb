#!/usr/bin/env ruby

require 'bio'	#	bioruby gem

if ARGV.empty?
	puts "\nUsage:"
	puts "\n#{$0} fasta_file\n\n"
	exit
end

filename  = ARGV[0]


#	while



inputfile     = Bio::FlatFile.auto(filename)
root_extname  = File.extname(filename)
root_filename = File.basename(filename,root_extname)

command = "grep '>' #{filename} | wc -l"
puts "Counting sequences with ..."
puts command
total_sequences = `#{command}`.to_i
#		#	will return string with leading spaces and trailing carriage return.  
#		#=> " 8899236\n"
#		#	to_i strips it all off and returns just the integer. Awesome.
#		##=> 8899236
#
puts "Found #{total_sequences} sequences"
#	using +1 as log of 2-10 is 1, 11-100 is 2
max_digits = Math.log10( total_sequences + 1 ).ceil

#names = {}
#File.open(name_file,'r').each{|name|
#	names[name.chomp] = 1
#}
#
#puts "Found #{names.length} sequence names."
#
puts "Opening output file."
#outputfile = File.open("#{name_file}.fasta",'w')
outputfile = File.open("#{root_filename}.uniq.fasta",'w')

sequences={}
duplicate_counter = 0
inputfile.each_with_index do |entry,index|
	printf "\rReading sequence %#{max_digits}d of %#{max_digits}d", index+1, total_sequences

#	if names.has_key?(entry.entry_id)
#		outputfile.puts entry
#	end

	if !sequences.has_key?(entry.seq)
		outputfile.puts entry
		sequences[entry.seq] = true
	else
		duplicate_counter += 1
	end

end
puts
puts "Removed #{duplicate_counter} duplicated sequences."
puts "Closing output files."
outputfile.close




#	end while

__END__
