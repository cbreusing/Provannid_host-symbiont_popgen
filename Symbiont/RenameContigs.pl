#!/usr/bin/env perl

use warnings;
use strict;

my $line;
my $index = 1;
my @outLines;

open (FILE, $ARGV[0]) or die "Cannot open input file: $!\n";

while (my $line = <FILE>)
 {
 if ($line =~ />/) {
 $line =~ s/>.*/>Scaffold_$index/g;
 $index = $index + 1;
 }
 push(@outLines,  $line);
 }
 
close FILE;
open (OUTFILE, ">$ARGV[1]");
print (OUTFILE @outLines);
close (OUTFILE);