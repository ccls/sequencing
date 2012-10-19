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

end
