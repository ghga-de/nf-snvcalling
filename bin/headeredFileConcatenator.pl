#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#


use strict;
use warnings;
use Pod::Usage;

my @files = @ARGV;
my $line;
my $colnames;
my $current_colnames;

if (!@ARGV) {
    pod2usage({ -verbose => 1, -message => "Error: No input files specified", -exitval => 2, -output => \*STDERR });
}

if ($ARGV[0] eq '--help') {
    pod2usage({ -verbose => 2, -exitval => 1 });
}

### First file: print header lines and set colnames
open(IN, $files[0]) || die "Could not open file $files[0] ($!)";

while ($line = <IN>) {
    print $line;
    last if ($line !~ /^\#/);
    $colnames = $line;
}
die "No header found in 1st file ($files[0])" if (!defined($colnames));

while (<IN>) {
    print;
}
close IN;

### Additional Files: do not print header; check if colnames match
foreach my $file (@files[1 .. $#files]) {
    open(IN, $file) || die "Could not open file $file ($!)";
    while ($line = <IN>) {
        last if ($line !~ /^\#/);
        $current_colnames = $line;
    }
    die "No header found in file $file" if (!defined($current_colnames));
    die "Columns in file $file do not match\n 1st file: $colnames\n Current file: $current_colnames\n" if ($colnames ne $current_colnames);
    print $line;

    while (<IN>) {
        print;
    }
    close IN;
}



__END__

=head1 NAME

headeredFileConcatenotor.pl

=head1 SYNOPSIS

headeredFileConcatenotor.pl file1 file2 ...

=head1 OPTIONS

=over 8

=item B<--help>

print help text

=back

=head1 DESCRIPTION

This script takes a number of files as command line arguments, concatenates them and print the concatenated file to stdout.
It prints out header lines (starting with #) from the 1st file but not from subsequent files.
It tests if the last header line (which should contain the column names) is identical in all files and exits with error if not.
Lines starting with # but are not at the beginning of the file (i.e. there were lines not starting with # in between) are not removed.

=cut
