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
open my $in, '<', $file;

my (@fields, $chr, $start, $end);

while (<$in>) {
    next if (/^\#/);
    @fields = split(/\t/);
    $chr = $fields[0];
    $start = ($fields[1] - $pad - 1 > 0) ? $fields[1] - $pad - 1 : 0;
    $end = $fields[1] + length($fields[3]) + $pad - 1; # TODO: prevent exceeding contig border
    say join "\t", $chr, $start, $end;
}
die "Error reading $file: $!\n" if $in->error;  # issue #43
close $in;
