#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#

# Fixed 2026-08-12 @kubranarci: Added autodie and guarded I/O reads with defined(readline())
# Changed behavior: Unguarded <FH> loops now detect read errors; I/O failures raise exceptions instead of silently returning undef

use strict;
use warnings;
use autodie;
use v5.10;

my $pad = shift;
my (@fields, $chr, $start, $end);
while (<>) {
    next if (/^\#/);
    @fields = split(/\t/);
    $chr = $fields[0];
    $start = ($fields[1] - $pad - 1 > 0) ? $fields[1] - $pad - 1 : 0;
    $end = $fields[1] + length($fields[3]) + $pad - 1; # TODO: prevent exceeding contig border
    say join "\t", $chr, $start, $end;
}
