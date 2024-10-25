#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 

use strict;
use warnings;
use v5.10;

# Retrieve output filenames from command-line arguments
my $snv_of = shift @ARGV or die "SNV output file not provided\n";
my $indel_of = shift @ARGV or die "Indel output file not provided\n";

# Open output files safely
open(my $snv_fh, '>', $snv_of) or die "Could not open $snv_of for writing: $!\n";
open(my $indel_fh, '>', $indel_of) or die "Could not open $indel_of for writing: $!\n";

# Process input from standard input
while (<>) {
    chomp;  # Remove newline at the end of the line
    if (/^\#/) {
        print $snv_fh "$_\n";  # Print headers to SNV file only
        next;
    }
    # Check if the 8th column indicates an indel
    if ((split(/\t/))[7] =~ /^INDEL/) {
        print $indel_fh "$_\n";
    } else {
        print $snv_fh "$_\n";
    }
}

# Close output files
close($snv_fh);
close($indel_fh);

