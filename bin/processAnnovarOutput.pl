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

(@ARGV >= 2) || die "Usage: processAnnovarOutput.pl <variant_function file> <exonic_variant_function file>";
my ($file1, $file2) = @ARGV;
my (@f2_line, @f1_line);
my $i;
my $nr = 0;
#  my $f2lnr; #file 2 line number
my @f2anno;
open(F1, $file1) || die "Could not open variant_function file";
open(F2, $file2) || die "Could not open exonic_variant_function file";
while (<F2>) {
    chomp;
    @f2_line = split(/\t/);
    #    $f2lnr = substr($line[0],4);
    @f2anno = ($f2_line[1], $f2_line[2]);
    while (<F1>) {
        chomp;
        #      $nr++;
        @f1_line = split(/\t/);
        if ($f2_line[3] eq $f1_line[2] && $f2_line[4] == $f1_line[3] && $f2_line[5] == $f1_line[4]) {
            say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], @f2anno;
            last;
        } else {
            say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
        }
    }
    die "Error reading variant_function file $file1: $!\n" if (\*F1)->error;  # issue #43
}
die "Error reading exonic_variant_function file $file2: $!\n" if (\*F2)->error;  # issue #43
while (<F1>) {
    # print out remaining lines from file 1 after file 2 has ended
    chomp;
    @f1_line = split(/\t/);
    say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
}
die "Error reading variant_function file $file1: $!\n" if (\*F1)->error;  # issue #43
close(F1);
close(F2);
