require 'test/unit'


class Test::Unit::TestCase

	#	Merging of rails' "test" in lib/active_support/testing/declarative.rb 
	#	and my own "test_with_verbosity" method.  dots are so not verbose.
	def self.test(name, &block)
		test_name = "test_#{name.gsub(/\s+/,'_')}".to_sym
		defined = instance_method(test_name) rescue false
		raise "#{test_name} is already defined in #{self}" if defined
		if block_given?
#			define_method(test_name, &block)
			define_method(test_name) do
#				print "\n#{self.class.name.gsub(/Test$/,'').titleize} #{name}: "
				print "\n#{self.class.name.gsub(/Test$/,'')} #{name}: "
				yield
			end
		else
			define_method(test_name) do
				flunk "No implementation provided for #{name}"
			end
		end
	end

end
__END__
