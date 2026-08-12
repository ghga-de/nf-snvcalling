#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Fixed 2026-08-12 @kubranarci: Added autodie and guarded I/O reads with defined(readline())
# Changed behavior: Unguarded <FH> loops now detect read errors; I/O failures raise exceptions instead of silently returning undef

use strict;
use warnings;
use autodie;

my $infile = $ARGV[0];

if($infile =~ /\.gz/){open(IN, "zcat $infile |") or die "Could not open the $infile to detect the confidence column\n";}
else{open(IN, "<$infile") or die "Could not open the $infile to detect the confidence column\n";}

my $line;
while (defined($_ = readline(IN))) {
	chomp;
	$line=$_;
	last if($_ =~ /^#CHROM\s/);
}
close IN;
my $i = 0;
my @line = split("\t", $line);
while($i <= @line)
{
	last if($line[$i] =~ /^CONFIDENCE$/);
	$i++;
}
print $i;