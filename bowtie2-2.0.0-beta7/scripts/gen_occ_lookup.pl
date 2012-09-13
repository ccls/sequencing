#!/usr/bin/env perl -w

#
# Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Generate lookup table that, given a packed DNA byte (four bases) and
# a character (A, C, G or T), returns how many times that character
# occurs in that packed byte.  Useful for quickly counting character
# occurrences in long strings.  The LUT is indexed first by character
# (0-3) then by byte (0-255).
#
# Larger lookup tables are also possible, though they seem
# counterproductive.  E.g., looking up eight bases at a time yields a
# 256K LUT, which doesn't fit in L1.  A four-base LUT is 1KB, easily
# fitting in L1.
#
# See ebwt.h. 
#

my @as4 = (), @as3 = (), @as2 = (), @as1 = ();
my @cs4 = (), @cs3 = (), @cs2 = (), @cs1 = ();
my @gs4 = (), @gs3 = (), @gs2 = (), @gs1 = ();
my @ts4 = (), @ts3 = (), @ts2 = (), @ts1 = ();

# Compile character arrays
my $i;
for($i = 0; $i < 256; $i++) {
	my $b01 = ($i >> 0) & 3;
	my $b23 = ($i >> 2) & 3;
	my $b45 = ($i >> 4) & 3;
	my $b67 = ($i >> 6) & 3;
	
	my $a4 = ($b01 == 0) + ($b23 == 0) + ($b45 == 0) + ($b67 == 0);
	my $c4 = ($b01 == 1) + ($b23 == 1) + ($b45 == 1) + ($b67 == 1);
	my $g4 = ($b01 == 2) + ($b23 == 2) + ($b45 == 2) + ($b67 == 2);
	my $t4 = ($b01 == 3) + ($b23 == 3) + ($b45 == 3) + ($b67 == 3);

	push @as4, $a4;
	push @cs4, $c4;
	push @gs4, $g4;
	push @ts4, $t4;

	my $a3 = ($b01 == 0) + ($b23 == 0) + ($b45 == 0);
	my $c3 = ($b01 == 1) + ($b23 == 1) + ($b45 == 1);
	my $g3 = ($b01 == 2) + ($b23 == 2) + ($b45 == 2);
	my $t3 = ($b01 == 3) + ($b23 == 3) + ($b45 == 3);

	push @as3, $a3;
	push @cs3, $c3;
	push @gs3, $g3;
	push @ts3, $t3;

	my $a2 = ($b01 == 0) + ($b23 == 0);
	my $c2 = ($b01 == 1) + ($b23 == 1);
	my $g2 = ($b01 == 2) + ($b23 == 2);
	my $t2 = ($b01 == 3) + ($b23 == 3);

	push @as2, $a2;
	push @cs2, $c2;
	push @gs2, $g2;
	push @ts2, $t2;

	my $a1 = ($b01 == 0) + 0;
	my $c1 = ($b01 == 1) + 0;
	my $g1 = ($b01 == 2) + 0;
	my $t1 = ($b01 == 3) + 0;

	push @as1, $a1;
	push @cs1, $c1;
	push @gs1, $g1;
	push @ts1, $t1;
}

my $entsPerLine = 16;

print "#include <stdint.h>\n\n";
print "/* Generated by gen_lookup_tables.pl */\n\n";

# Count occurrences in all 4 bit pairs 

print "uint8_t cCntLUT_4[4][4][256] = {\n";
print "\t/* All 4 bit pairs */ {\n";

# Print As array
print "\t\t/* As */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$as4[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Cs array
print "\t\t/* Cs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$cs4[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Gs array
print "\t\t/* Gs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$gs4[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Ts array
print "\t\t/* Ts */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$ts4[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t}\n\t},\n";

# Count occurrences in low 1 bit pair

print "\t/* Least significant 1 bit pair */ {\n";

# Print As array
print "\t\t/* As */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$as1[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Cs array
print "\t\t/* Cs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$cs1[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Gs array
print "\t\t/* Gs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$gs1[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Ts array
print "\t\t/* Ts */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$ts1[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t}\n\t},\n";

# Count occurrences in low 2 bit pairs

print "\t/* Least significant 2 bit pairs */ {\n";

# Print As array
print "\t\t/* As */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$as2[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Cs array
print "\t\t/* Cs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$cs2[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Gs array
print "\t\t/* Gs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$gs2[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Ts array
print "\t\t/* Ts */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$ts2[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t}\n\t},\n";

# Count occurrences in low 3 bit pairs

print "\t/* Least significant 3 bit pairs */ {\n";

# Print As array
print "\t\t/* As */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$as3[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Cs array
print "\t\t/* Cs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$cs3[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Gs array
print "\t\t/* Gs */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$gs3[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t},\n";

# Print Ts array
print "\t\t/* Ts */ {\n";
for($i = 0; $i < 256; $i++) {
	print "\t\t\t" if(($i % $entsPerLine) == 0);
	print "$ts3[$i], ";
	print "\n" if(($i % $entsPerLine) == ($entsPerLine-1));
}
print "\t\t}\n\t}\n";

print "};\n";
