#!/usr/bin/env perl -w

use strict;
use warnings FATAL => 'all';

#do_this('echo "begin";sleep 4; echo "end"; funkymonky');

#do_this('touch my_empty_file');

#print 'exist' if( -e 'missing_file' );
#print -s 'missing_file';
#print -s 'my_empty_file';
#print -s $0;



#check_file( 'my_empty_file' );

#my ($sec,$min,$hour,$mday,$mon,$year,
#          $wday,$yday,$isdst) = localtime time;
#my $dstring = sprintf( "%4d%02d%02d%02d%02d%02d.outdir",
#	$year+1900, $mon+1, $mday, $hour, $min, $sec );
#
#print $dstring;
#mkdir $dstring;
#chdir $dstring;


#my $data = ">comp1810_c0_seq1 len=430 path=[408:0-429]";
##my ($front, $back) = split (/len[:]/, $data);
##my ($front, $back) = split (/len[=]/, $data);
#my ($front, $back) = split (/len=/, $data);
#print "$front .... $back\n"; # there is no "len[:]" in $data[0] !!!

#if(( 1 > 0 ) && ( 2 < 3 )){
#	print "yes";
#}

sub check_file {
	my ( $filename, $size ) = @_;
	$size ||= 10;
print "checking $filename\n";
	die "$filename not created" unless -e $filename;
print "$filename exists.\n";
	die "$filename empty ( <= $size )" if( -s $filename <= $size );
print "$filename not empty.\n";
}



sub do_this {
	my ( $command ) = @_;
	print ( "Executing ...\n" );
	print ( "$command\n" );

#	backticks cache the output so you see nothing until its done.
#	@result would be that output that COULD be printed.
#	my @result = `$command`;	

#	$exit_status here is the same as the built in $?
#	so long as another system command hasn't been run.
#	my $exit_status = system($command);
	system($command);

	print ( "\n$command failed with $?\n\n" ) if ( $? );
} 


__END__

Table 11-1.  File tests and their meanings
File test 	Meaning
-r 	File or directory is readable by this (effective) user or group
-w 	File or directory is writable by this (effective) user or group
-x 	File or directory is executable by this (effective) user or group
-o 	File or directory is owned by this (effective) user
-R 	File or directory is readable by this real user or group
-W 	File or directory is writable by this real user or group
-X 	File or directory is executable by this real user or group
-O 	File or directory is owned by this real user
-e 	File or directory name exists
-z 	File exists and has zero size (always false for directories)
-s 	File or directory exists and has nonzero size (the value is the size in bytes)
-f 	Entry is a plain file
-d 	Entry is a directory
-l 	Entry is a symbolic link
-S 	Entry is a socket

File Test Operators

Table 11-1.  File tests and their meanings (continued)
File test 	Meaning
-p 	Entry is a named pipe (a “fifo”)
-b 	Entry is a block-special file (like a mountable disk)
-c 	Entry is a character-special file (like an I/O device)
-u 	File or directory is setuid
-g 	File or directory is setgid
-k 	File or directory has the sticky bit set
-t 	The filehandle is a TTY (as reported by theisatty()system function; filenames can’t be tested by this test)
-T 	File looks like a “text” file
-B 	File looks like a “binary” file
-M 	Modification age (measured in days)
-A 	Access age (measured in days)
-C 	Inode-modification age (measured in days) 
Read more at http://www.devshed.com/c/a/Perl/File-Tests-in-Perl/#PZVE8rlys2YLSC1q.99

