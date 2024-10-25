#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 


use strict;
use warnings;
use v5.10;

(@ARGV >= 2) || die "Usage: processAnnovarOutput.pl <variant_function file> <exonic_variant_function file>";
my ($file1, $file2) = @ARGV;
my (@f2_line, @f1_line);
my $i;
my $nr = 0;
#  my $f2lnr; #file 2 line number
my @f2anno;
open(my $f1_fh, '<', $file1) || die "Could not open variant_function file '$file1': $!";
open(my $f2_fh, '<', $file2) || die "Could not open exonic_variant_function file '$file2': $!";

while (my $line2 = <$f2_fh>) {
    chomp $line2;
    my @f2_line = split(/\t/, $line2);
    my @f2_anno = ($f2_line[1], $f2_line[2]);

    # Reset the file1 handle to the beginning for each line of file2
    seek($f1_fh, 0, 0); 

    my $match_found = 0;
    while (my $line1 = <$f1_fh>) {
        chomp $line1;
        my @f1_line = split(/\t/, $line1);

        if ($f2_line[3] eq $f1_line[2] && $f2_line[4] == $f1_line[3] && $f2_line[5] == $f1_line[4]) {
            say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], @f2_anno;
            $match_found = 1;
            last;  # Exit the inner loop upon finding a match
        }
    }

    # If no match was found, print the output with placeholders
    if (!$match_found) {
        say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
    }
}
# Print remaining lines from file1 after file2 has ended
while (my $line1 = <$f1_fh>) {
    chomp $line1;
    my @f1_line = split(/\t/, $line1);
    say join "\t", $f1_line[2], $f1_line[7], $f1_line[8], $f1_line[0], $f1_line[1], '.', '.';
}

close($f1_fh);
close($f2_fh);
