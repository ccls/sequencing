class String

	def execute
		puts "Executing ..."
		puts self
		#	status is true or false
		status = system(self);
		#	$? is the pid and the return code
		puts ( "\n#{self} failed with #{$?}\n" ) unless ( status );
	end

	#
	#	A leading @ would be from a FASTQ file.
	#	RINS' fastq2fasta.pl script DOES NOT remove the leading @
	#	This is NOT part of the sequence name and when data makes
	#		it into a SAM file, it will fail.
	#	I modified our version to remove them, but also do so here
	#		just in case.
	#
	#	A leading > would be from a FASTA file.
	#
	#	Trailing /1 or /2 is from Illumina
	#
	#	Suffixes like _2:N:0:CGTAGG are Casava 1.8 (modified)
	#	the fastq2fasta.pl script replaces a SPACE with the _
	#
	def delane_sequence_name
		name = self.chomp
		# >@HWI-ST281_0133:3:1:1222:2139#0/1
		name.gsub!(/^>/,'')
		# @HWI-ST281_0133:3:1:1222:2139#0/1
		# 
		# @HWI-ST281_0133:3:1:1222:2139#0/1
		name.gsub!(/^@/,'')
		# HWI-ST281_0133:3:1:1222:2139#0/1
		# 
		# HWI-ST281_0133:3:1:1222:2139#0/1
		name.gsub!(/\/\d+$/,'')
		# HWI-ST281_0133:3:1:1222:2139#0
		# 
		# HWI-ST977:132:C09W8ACXX:7:2307:10304:95858_2:N:0:CGTAGG
		name.gsub!(/_\d{1}:.*$/,'')
		# HWI-ST977:132:C09W8ACXX:7:2307:10304:95858
		# 
		name
	end


	def file_check( die_on_failed_file_check = false, empty_size = 0 )
		msg = '';
		unless( File.exists?(self) )
			msg = "#{self} not created";
			if( die_on_failed_file_check )
				raise msg;
			else
				puts msg
			end
		end
		#	Need to check existance too because if not dying on failure
		#	could have passed above check.  Then '-s $filename' is '' raising ...
		#	Use of uninitialized value in numeric le (<=) at ...
		if( ( File.exists?(self) ) && ( File.size(self) <= empty_size ) )
			msg = "#{self} empty ( <= #{empty_size} )";
			if( die_on_failed_file_check )
				raise msg
			else
				puts msg
			end
		end
	end

end

#
#	Called from blatoutcandidate.rb and pull_reads_fasta.rb
#
def pull_reads( names, input_fasta, output_fasta)
	File.open( input_fasta,'r') { |input|
	File.open(output_fasta,'w') { |output|
		while line = input.gets do
			if( names[line.delane_sequence_name] )
				output.puts line
				output.puts input.gets
			else
				input.gets
			end
		end
	} }
end

class CclsSequencer

	attr_accessor :options
	attr_accessor :original_stdout
	attr_accessor :original_stderr

	def initialize(options={})
		self.options = {
#			:config_filename          => 'config.yml',	#	used in the app, not here
			:output_filename          => 'results.txt',
			:output_suffix            => 'dark',
			:bowtie_version           => 1,				#	irrelevant now (dark uses 2, rins uses 1)
			:link_sample_fa_files     => false,		#	used in file_format_check_and_conversion (rins only)
			:pre_chopped              => false,		#	rins only
			:chop_read_length         => 25,			#	rins only
			:minIdentity              => 80,
			:compress_ratio_thrd      => 0.5,
			:iteration                => 2,				#	rins only
			:bowtie_threads           => 4,	#6,
			:bowtie_mismatch          => 3,				#	rins only (bowtie2 uses capital N which can only be 0 or 1)
			:paired_fragment_length   => 300,
			:min_contig_length        => 300,
			:trinity_threads          => 6,
			:blastn_evalue_thrd       => 0.05,
			:blastn_outfmt            => 6,
			:similarity_thrd          => 0.8,
			:mailto                   => '',			#	no longer used
			:die_on_failed_file_check => false,
			:file_format              => 'fasta',
			:blat_reference           => ['/Volumes/cube/working/indexes/all_bacterial.fna',
			                              '/Volumes/cube/working/indexes/virus_all.fasta' ],
			:bowtie_index_human       => '/Volumes/cube/working/indexes/hg18',
			:bowtie_human_indexes     => ['/Volumes/cube/working/indexes/hg18',
			                              '/Volumes/cube/working/indexes/hg19',
			                              '/Volumes/cube/working/indexes/Blast1',
			                              '/Volumes/cube/working/indexes/Blast2',
			                              '/Volumes/cube/working/indexes/Homo_sapiens.GRCh37.69.cdna.all'],
			:blastn_index_human       => '/Volumes/cube/working/indexes/hg18',
			:blastn_index_non_human   => '/Volumes/cube/working/indexes/nt',
			:files                    => {}
		}.merge(options)
	end

	def method_missing(symb,*args,&block)
		#
		#	basically, check the config options for key/value
		#	so don't need to 'options[key]' and can just 'key'
		#
		if options.has_key? symb.to_s.to_sym
			options[symb.to_s.to_sym]
		else
			super
		end
	end

	def prepare_output_dir_and_log_file
		puts "Preparing working dir and log file"
		outdir = ["#{Time.now.strftime("%Y%m%d%H%M%S")}.outdir",
			options[:output_suffix]].compact.join('.')
		puts "Working dir is #{outdir}"
		FileUtils.mkdir outdir
		FileUtils.chdir outdir

		puts "About to redirect STDOUT and STDERR."
		puts "No output should go to the screen until complete."
		puts "Perhaps put this process in the background with Ctrl-Z then 'bg'."
		puts "Then use the following command to follow along ..."
		puts "tail -f #{outdir}/log_file.txt"
		self.original_stdout = STDOUT.clone
		self.original_stderr = STDERR.clone
		STDOUT.reopen('log_file.txt','a')
		STDERR.reopen('log_file.txt','a')
	end

	def file_format_check_and_conversion
#
#	TODO use my file_format_detector instead of the value in the config file?
#
#		file_format = FileFormatDetector.new(files.first.value).format
#

		raise "File format can either be fastq or fasta" unless( 
			['fasta','fastq'].include?(file_format) )
		
		puts "step 1 change fastq files to fasta files"
		files.each_pair do |k,v|
			if( file_format == "fastq")
				"fastq2fasta.pl #{v} #{k}lane.fa".execute
				"#{k}lane.fa".file_check(die_on_failed_file_check)
			else #if ($file_format eq "fasta") {
				if( link_sample_fa_files )
					puts "already fasta format, linking #{v} #{k}lane.fa instead"
					FileUtils.ln_s(v,"#{k}lane.fa")
					"#{k}lane.fa".file_check(die_on_failed_file_check)
				else
					puts "already fasta format, copying #{v} #{k}lane.fa instead"
					FileUtils.cp(v,"#{k}lane.fa")
					"#{k}lane.fa".file_check(die_on_failed_file_check)
				end
			end
		end
	end

	def pull_reads_from_fastas(names, fastas, outputs)
		puts "pull reads from #{fastas.join(',')} fasta files"
		puts "(join and unique laneless sequence names from the names files"
		puts " and select those from the inputs and place in the outputs.)"
		command = "pull_reads_fasta.rb "
		command << "#{names.join(' ')} "
		command << "#{fastas.join(' ')} "
		command << "#{outputs.join(' ')} "
		command.execute
		outputs.each{|output| output.file_check(die_on_failed_file_check) }
	end
#
#	These two and VERY similar.
#
	def blat_out_candidate_reads(blatpsls, fastas, outputs)
		puts "find blat out candidate reads"
		puts "(join and unique laneless sequence names from psl files"
		puts " and select those from the inputs and place in the outputs.)"
		command = "blatoutcandidate.rb "
		command << "#{blatpsls.join(' ')} "
		command << "#{fastas.join(' ')} "
		command << "#{outputs.join(' ')} "
		command.execute
		outputs.each{|output| output.file_check(die_on_failed_file_check) }
	end

	def blastn_non_human(fasta)
		rootname = fasta.gsub(/#{File.extname(fasta)}$/,'')
		command = "blastn -query=#{fasta} -db=#{blastn_index_non_human} "
		command << "-evalue #{blastn_evalue_thrd} " #unless blastn_evalue_thrd.blank?
		command << "-outfmt #{blastn_outfmt} " #unless blastn_outfmt.blank?
		command << "> #{rootname}_blastn.txt"
		command.execute
		"#{rootname}_blastn.txt".file_check(die_on_failed_file_check)
	end
		
	def wrap_things_up
		puts "Finished at ..."
		system("date")
		STDOUT.reopen(original_stdout)
		STDERR.reopen(original_stderr)
		puts "All done."
		system("date")
	end

end
