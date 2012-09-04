#!/usr/bin/env rake

require 'rake'
require 'rake/testtask'
require 'rdoc/task'

desc 'Default: run unit tests.'
task :default => :test

desc 'Test the html_test plugin.'
Rake::TestTask.new(:test) do |t|
        t.libs << 'test'
        t.libs << 'lib'
        t.pattern = 'test/**/*_test.rb'
        t.verbose = true
end
