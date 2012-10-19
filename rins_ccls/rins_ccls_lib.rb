class String

	def execute
		puts "Executing ..."
		puts self
		#	status is true or false
		status = system(self);
		#	$? is the pid and the return code
		puts ( "\n#{self} failed with #{$?}\n" ) unless ( status );
	end

	def delane_sequence_name
		name = self.chomp
		name.gsub!(/^>/,'')
		# 
		# >@HWI-ST281_0133:3:1:1222:2139#0/1
		name.gsub!(/\/\d+$/,'')
		# >@HWI-ST281_0133:3:1:1222:2139#0
		# 
		# > @HWI-ST977:132:C09W8ACXX:7:2307:10304:95858_2:N:0:CGTAGG
		name.gsub!(/_\d{1}:.*$/,'')
		# > @HWI-ST977:132:C09W8ACXX:7:2307:10304:95858
		# 
		name
	end


	def file_check( die_on_failed_file_check = false, empty_size = 0 )
		filename = self
	
	#	STDERR is not logged when using " | tee -a log"
	
		msg = '';
		unless( File.exists?(filename) )
			msg = "#{filename} not created";
			if( die_on_failed_file_check )
				raise msg;
			else
				puts msg
			end
		end
		#	Need to check existance too because if not dying on failure
		#	could have passed above check.  Then '-s $filename' is '' raising ...
		#	Use of uninitialized value in numeric le (<=) at ...
		if( ( File.exists?(filename) ) && ( File.size(filename) <= empty_size ) )
			msg = "#{filename} empty ( <= #{empty_size} )";
			if( die_on_failed_file_check )
				raise msg
			else
				puts msg
			end
		end
	end

end

class CclsSequencer

	attr_accessor :options
	attr_accessor :original_stdout
	attr_accessor :original_stderr

	def initialize(options={})
		self.options = options
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

	def wrap_things_up
		puts "Finished at ..."
		system("date")
		STDOUT.reopen(original_stdout)
		STDERR.reopen(original_stderr)
		puts "All done."
		system("date")
	end

end
