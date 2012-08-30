#!/usr/bin/env ruby

require 'csv'
require 'fileutils'

csv_file = 'Regions to verify.csv'

(f=CSV.open( csv_file, 'rb',{ :headers => true })).each do |line|
	puts
	puts "Processing line :#{f.lineno}:"
	puts line

	bed_files = Dir["#{line['SAMPLE ID']}*bed"]
	raise "No files matching :#{line['SAMPLE ID']}*bed: found" if bed_files.length < 1
	raise "Multiple files matching :#{line['SAMPLE ID']}*bed: found" if bed_files.length > 1
	
	bed_files.each do |bed_file_name|
		puts '-' * 70
		out_file_basename = "#{line['SAMPLE ID']}_#{line['#Chr1']}_#{line['Pos1']}"
		puts "Opening #{out_file_basename} fastq files for writing"
		FileUtils.mkdir("out") unless Dir.exists?("out")
		File.open("out/#{out_file_basename}.1.fastq",'w') do |out_1_file|
		File.open("out/#{out_file_basename}.2.fastq",'w') do |out_2_file|
			puts "Reading :#{bed_file_name}"
			File.open(bed_file_name,'r') do |bed_file|
				flag = false
				bed_line_no = 0
				while bed_line = bed_file.gets
					bed_line_no += 1
#track name=chr11_92083328_INV_673	description="BreakDancer chr11 92083328 INV 673"	useScore=0
					puts "#{bed_line_no}:#{bed_line}"
					if bed_line.match(/track name=/)
						puts "Found track line"
						if( bed_line.match(/track name=#{line['#Chr1']}_#{line['Pos1']}/) ) 
							puts "---- Matches :#{line['#Chr1']}_#{line['Pos1']}: setting flag to TRUE"
							flag = true 
						else
							puts "no match :#{line['#Chr1']}_#{line['Pos1']}: setting flag to false"
							flag = false
						end
					elsif flag
						puts "flag is true.  parsing line for HWI and sign"

#find HWI-ST1215:67:C0R7FACXX:4:2203:18517:186061 in ... fastq files (1 or 2) should only be one (but processing both)


# all of the output is duplicated because all of the HWI numbers are duplicated
#chrchr9	37090827	37090927	HWI-ST1215:67:C0R7FACXX:6:2106:11132:18849|JoeWiem_hg19_req	370	-	37090827	37090927	0,0,255
#chrchr9	37090817	37090917	HWI-ST1215:67:C0R7FACXX:6:2106:11132:18849|JoeWiem_hg19_req	370	+	37090817	37090917	255,0,0
#	one with a + and one with a -
#	as I am ignoring them, I search for and find them twice.

#	I can either be more specific and look for the + or - 
#	and then search the appropriate file
#	or just remember that I already searched for one.


#	+  corresponds to the 1.fastq file
#	-  corresponds to the 2.fastq file
	
						bed_line_array = bed_line.match(/(HWI-.*)\|.*(\+|\-)/)
						hwi = bed_line_array[1]
						sign = bed_line_array[2]
						raise "No HWI found in:#{bed_line}:" unless hwi
						raise "No sign found in:#{bed_line}:" unless sign

						file_name_num, out_file = case sign 
							when '+' then [1, out_1_file]
							when '-' then [2, out_2_file]
							else raise "Invalid sign :#{sign}:"
						end
						puts "Found HWI:#{hwi}: and sign:#{sign}:"	

						#	MAKE SURE TO NOT READ THIS FILES GENERATED FASTQ FILES! (added .remdup)
						#	Doesn't matter much now as I am writing to the 'out' directory
						fastq_files = Dir["#{line['SAMPLE ID']}*#{file_name_num}.fastq"]
						raise "No fastq files matching :#{line['SAMPLE ID']}*#{file_name_num}.fastq: found" if fastq_files.length < 1
						raise "Multiple fastq files matching :#{line['SAMPLE ID']}*#{file_name_num}.fastq: found" if fastq_files.length > 1

						fastq_files.each do |fastq_file_name|
							puts "Reading fastq file:#{fastq_file_name}:"
							fastq_file_contents = IO.read(fastq_file_name)
#							out_file.puts fastq_file_contents.match(/^@#{hwi}\n\w+\n\+\n.+\n/)
#	match finds first, scan finds all (but should only be one match anyway)
							matches = fastq_file_contents.scan(/^@#{hwi}\n\w+\n\+\n.+\n/)
							puts "Found #{matches.length}"
raise "No matches???" if matches.length < 1
							out_file.puts matches
#							puts fastq_file_contents.scan(/^@#{hwi}\n\w+\n\+\n.+\n/)
#							matches = fastq_file_contents.match(/^@#{hwi}\n\w+\n\+\n.+\n/)
						end	#	Dir["#{line['SAMPLE ID']}*fastq"].each do |fastq_file_name|
#memory["#{line['SAMPLE ID']}#{hwi}"] = true



					end	#	elsif flag

				end	#	while line = bed_file.gets
			end	#	File.open(bed_file_name) do |bed_file|
		end	#	File.open("#{out_file_basename}.1.fastq",'w') do |out_1_file|
		end	#	File.open("#{out_file_basename}.2.fastq",'w') do |out_2_file|
	end	#	Dir["#{line['SAMPLE ID']}*bed"].each do |bed_file|
end	#	(f=CSV.open( csv_file, 'rb',{ :headers => true })).each do |line|
