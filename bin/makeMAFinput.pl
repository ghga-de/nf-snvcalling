#!/usr/bin/env perl

# tumor variant frequency for high confidence somatic SNVs as input for MAF plots
# output: PID (number!)	chrom	tumor variant frequency
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 

use strict;
use warnings;

if (@ARGV < 2) {
    die "1. VCF file to calculate tumor variant frequency for MAF plot\n" .
        "2. Minimal coverage cutoff (default 0)\n" .
        "3. Minimal confidence threshold (default 8)\n";
}

# Get the input file
my $file = shift;

# Open the filehandle safely
open(my $fh, '<', $file) or die "Could not open file '$file': $!\n";

# Initialize header variable and found flag
my $header = "";
my $found = 0;

# Read through the file to find the header line
while ($header = <$fh>) {
    # Exit the loop if we find the '#CHROM' line (the line with the column names)
    if ($header =~ /^\#CHROM/) {
        $found = 1;  # Mark that we found the header.
        last;
    }
}

# Handle case where header was not found
if (!$found) {
    die "Header line not found in file: $file\n";
}
chomp $header;

my $CONFIDENCE = 28;
my $INFO = 7;
my $DBSBP = 0;

my @help = split(/\t/, $header);
for (my $c = 0; $c < @help; $c++)
{
	# fixed column labels (from vcf_pileup_compare)
	if ($help[$c] eq "CONFIDENCE")
	{
		$CONFIDENCE = $c;
		print STDERR "CONFIDENCE in column $c\n";
	}
	if ($help[$c] eq "INFO")
	{
		$INFO = $c;
		print STDERR "INFO in column $c\n";
	}
	if ($help[$c] eq "DBSNP")
	{
		$DBSBP = $c;
		print STDERR "DBSNP_COL in column $c\n";
	}
}

###############################################################################


my $mincov = shift;
unless (defined $mincov)
{
	$mincov = 0;
}
my $minconfidence = shift;
unless (defined $minconfidence)
{
	$minconfidence = 8;
}

#INFO => DP4 of tumor => depth and variant frequency
my ($trf, $trr, $tvf, $tvr) = (0, 0, 0, 0);
@help = ();
my $depth = 0;
# tumor variant frequency from DP4 fields (pos. 7 of vcf file)
my $tumvarfrq = 0;

# Assume that $CONFIDENCE, $INFO, $minconfidence, $mincov, $DBSBP are already defined earlier in the script.

# Loop through each line of the file
while (my $line = <$fh>)
{
    # Skip lines that start with '#'
    next if $line =~ /^#/;
    
    # Remove newline character at the end of the line
    chomp $line;
    
    # Split the line into fields based on tab character
    my @help = split("\t", $line);
    
    # Check if the confidence value meets the minimum threshold
    if ($help[$CONFIDENCE] >= $minconfidence) 
    {
        # Capture the four values (trf, trr, tvf, tvr) from the INFO field using a regex match
        if ($help[$INFO] =~ /DP4=(\d+),(\d+),(\d+),(\d+);/) 
        {
            my ($trf, $trr, $tvf, $tvr) = ($1, $2, $3, $4);
            
            # Calculate depth
            my $depth = $tvf + $tvr + $trf + $trr;
            
            # Check if depth meets the minimum coverage
            if ($depth >= $mincov) 
            {
                # Calculate the tumor variant frequency
                my $tumvarfrq = sprintf("%.2f", ($tvf + $tvr) / ($tvf + $tvr + $trf + $trr));
                
                # Check if there's an exact match in the DBSBP field and print the result
                if ($help[$DBSBP] =~ /MATCH=exact/) 
                {
                    print "1\t$help[0]\t$tumvarfrq\t$depth\n";
                } 
                else 
                {
                    print "0\t$help[0]\t$tumvarfrq\t$depth\n";
                }
            }
        } 
        else 
        {
            # Handle the case where the regex does not match
            warn "Line $.: INFO field does not match expected pattern: $line\n";
        }
    }
}

# Close the filehandle
close($fh) or warn "Could not close filehandle: $!";

exit;

