#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Fixed 2026-08-12 @kubranarci: I/O error detection (issue #43): autodie for open/close, and an explicit
# read-error check after each read loop (<FH> returns undef on a read error just like at EOF)

use strict;
use warnings;
use autodie;
use IO::Handle;

my $infile     = $ARGV[0];
my $from_stdin = (!defined $infile || $infile eq '-');
my $src        = $from_stdin ? 'STDIN' : $infile;

my $in_fh;
if ($from_stdin) {
    $in_fh = *STDIN;
} elsif ($infile =~ /\.gz$/) {
    open($in_fh, "zcat $infile |") or die "Could not open the $infile to detect the confidence column\n";
} else {
    open($in_fh, '<', $infile) or die "Could not open the $infile to detect the confidence column\n";
}

my $line;
my $found = 0;
while (<$in_fh>) {
    chomp;
    $line = $_;
    if (/^#CHROM\s/) { $found = 1; last }
}
die "Error reading $src: $!\n" if $in_fh->error;  # issue #43

unless ($from_stdin) {
    no autodie 'close';
    close $in_fh or warn "Could not close input handle for $src: $?";
}

die "No #CHROM header found in $src\n" unless $found;

my @line = split("\t", $line);
my ($i) = grep { $line[$_] eq 'CONFIDENCE' } 0 .. $#line;
die "No CONFIDENCE column found in $src\n" unless defined $i;
print $i;
