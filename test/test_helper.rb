require 'test/unit'


class String

	def namerize
#		self.downcase.gsub(/\b('?[a-z])/) { $1.capitalize }
#	what about O'Grady
		self.downcase.gsub(/\b([a-z])/) { $1.capitalize }
	end



# File activesupport/lib/active_support/inflector/methods.rb, line 94
#i	def humanize(lower_case_and_underscored_word)
	def humanize
		lower_case_and_underscored_word = self
		result = lower_case_and_underscored_word.to_s.dup
#		inflections.humans.each { |(rule, replacement)| break if result.gsub!(rule, replacement) }
		result.gsub!(/_id$/, "")
		result.gsub!(/_/, ' ')
#		result.gsub(/([a-z\d]*)/) { |match|
#		"#{inflections.acronyms[match] || match.downcase}"
#		}.gsub(/^\w/) { $&.upcase }
	end

# File activesupport/lib/active_support/inflector/methods.rb, line 77
#	def underscore(camel_cased_word)
	def underscore
		camel_cased_word = self
		word = camel_cased_word.to_s.dup
		word.gsub!(/::/, '/')
#		word.gsub!(/(?:([A-Za-z\d])|^)(#{inflections.acronym_regex})(?=\b|[^a-z])/) { "#{$1}#{$1 && '_'}#{$2.downcase}" }
		word.gsub!(/([A-Z\d]+)([A-Z][a-z])/,'\1_\2')
		word.gsub!(/([a-z\d])([A-Z])/,'\1_\2')
		word.tr!("-", "_")
		word.downcase!
		word
	end

# File activesupport/lib/active_support/inflector/methods.rb, line 115
#    def titleize(word)
	def titleize
#		humanize(underscore(self)).gsub(/\b('?[a-z])/) { $1.capitalize }
#		humanize(self.underscore).gsub(/\b('?[a-z])/) { $1.capitalize }
		self.underscore.humanize.gsub(/\b('?[a-z])/) { $1.capitalize }
	end

end

class Test::Unit::TestCase

	#	Merging of rails' "test" in lib/active_support/testing/declarative.rb 
	#	and my own "test_with_verbosity" method.  dots are so not verbose.
	def self.test(name, &block)
		test_name = "test_#{name.gsub(/\s+/,'_')}".to_sym
		defined = instance_method(test_name) rescue false
		raise "#{test_name} is already defined in #{self}" if defined
		if block_given?
			define_method(test_name,&block)
#
#	below works, but the tests can't print to the screen
#
#			define_method(test_name) do
#				print "\n#{self.class.name.gsub(/Test$/,'')} #{name}: "
#				#	need to wrap 'yield' in a lambda to delay calling
#				#	otherwise we get ....
##	NoMethodError: undefined method `assert' for UniquifySequenceNamesTest:Class
#				lambda{ yield }
#			end
		else
			define_method(test_name) do
				flunk "No implementation provided for #{name}"
			end
		end
	end



	def self.test_with_verbosity(name,&block)
		test_without_verbosity(name,&block)

		test_name = "test_#{name.gsub(/\s+/,'_')}".to_sym
		define_method("_#{test_name}_with_verbosity") do
			print "\n#{self.class.name.gsub(/Test$/,'').titleize} #{name}: "
#			print "\n#{self.class.name.gsub(/Test$/,'')} #{name}: "
			send("_#{test_name}_without_verbosity")
		end
		#
		#	can't do this...
		#		alias_method_chain test_name, :verbosity
		#	end up with 2 methods that begin
		#	with 'test_' so they both get run
		#
		alias_method "_#{test_name}_without_verbosity".to_sym,
			test_name
		alias_method test_name,
			"_#{test_name}_with_verbosity".to_sym
	end

	class << self
		alias_method :test_without_verbosity, :test
		alias_method :test, :test_with_verbosity
	end


end
