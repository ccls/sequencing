require 'test_helper'

class UniquifySequenceNamesTest < Test::Unit::TestCase

	test "should find duplicate sequence name" do
		outfile = 'test/data/duplicate_sequence_names.fa.out'	
		File.delete(outfile) if File.exists?(outfile)
		assert !File.exists?(outfile)
		system 'uniquify_sequence_names.rb test/data/duplicate_sequence_names.fa > /dev/null'

#puts
#puts IO.read(outfile)

		assert File.exists?(outfile)
		File.delete(outfile)
	end

end
