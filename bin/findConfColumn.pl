#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Modified 2026-08-12 @kubranarci: support '-' for STDIN, fail if no #CHROM header, fix off-by-one in column search

use strict;
use warnings;

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
    if ($_ =~ /^#CHROM\s/) {
        $line = $_;
        last;
    }
}

if ($infile && $infile ne '-') {
    # zcat gets SIGPIPE when we stop after the header; only report other failures
    close $in_fh or (($? & 127) == 13) or warn "Could not close input handle for $infile: exit status $?\n";
}

die "No #CHROM header found in $infile\n" unless defined $line;
my $i = 0;
my @line = split("\t", $line);
while($i <= $#line) {
    last if($line[$i] =~ /^CONFIDENCE$/);
    $i++;
}
print $i;