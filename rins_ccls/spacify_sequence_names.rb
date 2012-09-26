#!/usr/bin/env ruby

require 'optparse'

# This hash will hold all of the options
# parsed from the command-line by
# OptionParser.
options = {
	:dryrun => false
}

optparse = OptionParser.new do |opts|
	# Set a banner, displayed at the top of the help screen.
	opts.banner = "\nUsage: #{File.basename($0)} [options] " <<
			"fa_in_file1 fa_in_file2 ...\n\n" <<
		"The purpose of this script is to convert ....\n" <<
		">gi_294756_gb_L11599.1_ROTVP4X_Lamb_rotavirus_Lp14_(VP4)_gene\n" <<
		" ... into ...\n" <<
		">gi|294756|gb|L11599.1| ROTVP4X Lamb rotavirus Lp14 (VP4) gene\n" <<
		" ... replacing the first 3 _'s with |'s \n" <<
		" ... the 4th _ with '| '\n" <<
		" ... and the rest with spaces.\n\n"

	# Define the options, and what they do

	opts.on( '-d', '--dryrun', "Don't create output file" ) do
		options[:dryrun] = true
	end

	# This displays the help screen, all programs are assumed to have this option.
	#	Add extra "\n" to last option for aesthetics.
	opts.on( '-h', '--help', 'Display this screen',"\n" ) do
		puts opts
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
if ARGV.empty?
	puts optparse	#	Basically display the command line help
	exit
end

ARGV.each do |infilename|
	if File.exists? infilename
		puts "Processing #{infilename}"
	else
		puts "#{infilename} not found. Skipping."
		next
	end
	outfilename = "#{infilename}.#{File.basename($0)}.out"

	File.open(  infilename, 'r' ) { |infile|
	File.open( outfilename, 'w' ) { |outfile|
		while( line = infile.gets )
			if line.match(/^>/)
				puts "Found sequence name."
				puts line
#				line.gsub!(/\s*$/,'')
#				line.gsub!(/(_)[^_]+__/,'|')
#				line.gsub!(/_/,'|')

#				parts = line.split('__')
#				parts[0].gsub!(/_/,'|')
#				parts[1].gsub!(/_/,' ')
#				line = parts.join('| ')

#	ONE line in virus is different
#	>gi|320202719|ref|NC_014967.1| Honeysuckle ringspot virus, complete genome
#	need to differentiate from the likes of ...
# >gi_47177085_ref_NC_005874.1__Zinnia_leaf_curl_virus-assoc

				unless line.match(/\|/)
					parts = line.split(/_+/)
					line = parts[0..3].join('|')
					line << "| "
					line << parts[4..-1].join(' ')
				end

#				puts "Duplicate sequence name. Renaming."
				puts line
			end
#			puts line if options[:verbose]
			outfile.puts line unless options[:dryrun]
		end
	} }	#	File.open
end



__END__

jakewendt@fxdgroup-169-229-196-225 : sequencing 505> grep "^>" /Volumes/cube/working/indexes/virus.fa | head
>gi_220961751_gb_EU984101.1__Human_rotavirus_A_strain_0613158-CA_NSP4_(NSP4)_gene,_complete_cds
>gi_323513964_gb_HQ641354.1__Vibrio_phage_ICP1_2004_A,_complete_genome
>gi_323513733_gb_HQ641353.1__Vibrio_phage_ICP1_2001_A,_complete_genome
>gi_323513502_gb_HQ641352.1__Vibrio_phage_ICP1_2005_A,_complete_genome
>gi_323513275_gb_HQ641351.1__Vibrio_phage_ICP1_2006_A,_complete_genome
>gi_323513048_gb_HQ641350.1__Vibrio_phage_ICP1_2006_B,_complete_genome
>gi_323512820_gb_HQ641349.1__Vibrio_phage_ICP1_2006_C,_complete_genome
>gi_323512592_gb_HQ641348.1__Vibrio_phage_ICP1_2006_D,_complete_genome
>gi_323512298_gb_HQ641346.1__Vibrio_phage_ICP2_2006_A,_complete_genome
>gi_323512177_gb_HQ641344.1__Vibrio_phage_ICP3_2007_A,_complete_genome
jakewendt@fxdgroup-169-229-196-225 : sequencing 506> grep "^>" /Volumes/cube/working/indexes/all| head   
all_bacterial.00.nhr                 all_bacterial_and_virus.00.nhr 
all_bacterial.00.nin                 all_bacterial_and_virus.00.nin 
all_bacterial.00.nsq                 all_bacterial_and_virus.00.nsq 
all_bacterial.01.nhr                 all_bacterial_and_virus.01.nhr 
all_bacterial.01.nin                 all_bacterial_and_virus.01.nin 
all_bacterial.01.nsq                 all_bacterial_and_virus.01.nsq 
all_bacterial.fna                    all_bacterial_and_virus.fa 
all_bacterial.fna.original           all_bacterial_and_virus.fa.original 
all_bacterial.nal                    all_bacterial_and_virus.nal 
jakewendt@fxdgroup-169-229-196-225 : sequencing 506> grep "^>" /Volumes/cube/working/indexes/all_bacterial.fn| head
all_bacterial.fna           all_bacterial.fna.original  
jakewendt@fxdgroup-169-229-196-225 : sequencing 506> grep "^>" /Volumes/cube/working/indexes/all_bacterial.fna | head
>gi|158333233|ref|NC_009925.1| Acaryochloris marina MBIC11017, complete genome
>gi|158341503|ref|NC_009933.1| Acaryochloris marina MBIC11017 plasmid pREB8, complete sequence
>gi|158341621|ref|NC_009934.1| Acaryochloris marina MBIC11017 plasmid pREB9, complete sequence
>gi|158341329|ref|NC_009932.1| Acaryochloris marina MBIC11017 plasmid pREB7, complete sequence
>gi|158341140|ref|NC_009931.1| Acaryochloris marina MBIC11017 plasmid pREB6, complete sequence
>gi|158339871|ref|NC_009927.1| Acaryochloris marina MBIC11017 plasmid pREB2, complete sequence
>gi|158340643|ref|NC_009929.1| Acaryochloris marina MBIC11017 plasmid pREB4, complete sequence
>gi|158340917|ref|NC_009930.1| Acaryochloris marina MBIC11017 plasmid pREB5, complete sequence
>gi|158340280|ref|NC_009928.1| Acaryochloris marina MBIC11017 plasmid pREB3, complete sequence
>gi|158339488|ref|NC_009926.1| Acaryochloris marina MBIC11017 plasmid pREB1, complete sequence
