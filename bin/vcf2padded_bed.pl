#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#

# Fixed 2026-08-12 @kubranarci: I/O error detection (issue #43): autodie for open/close, and an explicit
# read-error check after each read loop (<FH> returns undef on a read error just like at EOF)

use strict;
use warnings;
use autodie;
use IO::Handle;
use v5.10;

my $pad  = shift;
my $file = shift;
my $fai  = shift;   # optional: reference .fai, to keep intervals inside the contig (#45)
open my $in, '<', $file;

my %contig_length;
if (defined $fai) {
    open my $fai_fh, '<', $fai;
    while (<$fai_fh>) {
        my ($name, $length) = split(/\t/);
        $contig_length{$name} = $length;
    }
    die "Error reading $fai: $!\n" if $fai_fh->error;
    close $fai_fh;
}

my (@fields, $chr, $start, $end);

while (<$in>) {
    next if (/^\#/);
    @fields = split(/\t/);
    $chr = $fields[0];
    $start = ($fields[1] - $pad - 1 > 0) ? $fields[1] - $pad - 1 : 0;
    $end = $fields[1] + length($fields[3]) + $pad - 1;
    # bedtools skips intervals that go beyond the contig end (#45)
    $end = $contig_length{$chr} if (exists $contig_length{$chr} && $end > $contig_length{$chr});
    say join "\t", $chr, $start, $end;
}
die "Error reading $file: $!\n" if $in->error;  # issue #43
close $in;
