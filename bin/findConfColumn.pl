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
my $in_fh;

if (!defined $infile || $infile eq '-') {
    $in_fh = *STDIN;
} elsif ($infile =~ /\.gz$/) {
    open($in_fh, "zcat $infile |") or die "Could not open the $infile to detect the confidence column\n";
} else {
    open($in_fh, '<', $infile) or die "Could not open the $infile to detect the confidence column\n";
}

my $line;
while (defined($_ = readline($in_fh))) {
    chomp;
    $line=$_;
    last if($_ =~ /^#CHROM\s/);
}

if ($infile && $infile ne '-') {
    no autodie 'close';
    close $in_fh or warn "Could not close input handle for $infile: $?";
}

die "No #CHROM header found in $infile\n" unless defined $line;
my $i = 0;
my @line = split("\t", $line);
while($i <= $#line) {
    last if($line[$i] =~ /^CONFIDENCE$/);
    $i++;
}
print $i;