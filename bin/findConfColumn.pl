#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 

use strict;
use warnings;

my $infile = shift @ARGV or die "Usage: $0 <input_file>\n";
die "File does not exist: $infile\n" unless -e $infile;
die "File is not readable: $infile\n" unless -r $infile;

# Open the input file safely
my $in_fh;  # Declare a lexical filehandle
if ($infile =~ /\.gz$/) {
    open($in_fh, "zcat $infile |") or die "Could not open gzipped file $infile: $!\n";
} else {
    open($in_fh, '<', $infile) or die "Could not open file $infile: $!\n";
}

my $line;
while (<$in_fh>) {
    chomp;
    $line = $_;
    last if $line =~ /^#CHROM\s/;
}

close $in_fh;  # Close the filehandle when done

my $i = 0;
my @line = split("\t", $line);
while($i <= @line)
{
	last if($line[$i] =~ /^CONFIDENCE$/);
	$i++;
}
print $i;