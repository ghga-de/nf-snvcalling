#!/usr/bin/env perl
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 


use strict;
use warnings;
use Pod::Usage;

# Check for input files
my @files = @ARGV;

if (!@files) {
    pod2usage({ -verbose => 1, -message => "Error: No input files specified", -exitval => 2, -output => \*STDERR });
}

if ($files[0] eq '--help') {
    pod2usage({ -verbose => 2, -exitval => 1 });
}

# Function to process a file and return the header
sub process_file {
    my ($filename, $print_header) = @_;
    my $current_colnames;
    
    open(my $fh, '<', $filename) or die "Could not open file $filename: $!\n";

    # Read lines and capture header
    while (my $line = <$fh>) {
        chomp($line);
        if ($line =~ /^\#/) {
            print "$line\n" if $print_header;  # Print header if requested
            $current_colnames = $line;          # Store header
            next;
        }
        last;  # Exit loop when a non-header line is encountered
    }
    
    die "No header found in file $filename\n" if !defined($current_colnames);
    
    return ($current_colnames, $fh);  # Return the header and the filehandle
}

# Process the first file
my ($colnames, $fh_first) = process_file($files[0], 1);  # Print header

# Read the rest of the files and check headers
foreach my $file (@files[1..$#files]) {
    my ($current_colnames, $fh_current) = process_file($file, 0);  # Do not print header
    
    # Compare headers
    die "Columns in file $file do not match\nFirst file: $colnames\nCurrent file: $current_colnames\n"
        if $colnames ne $current_colnames;

    # Print the rest of the current file
    while (my $line = <$fh_current>) {
        print $line;
    }
    
    close $fh_current;  # Close the current filehandle
}

close $fh_first;  # Close the first filehandle





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
