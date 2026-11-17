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

my $fastaBinary = $ARGV[0];
my $file = $ARGV[1];
my $refGenome = $ARGV[2];
my $pad = $ARGV[3];
my %fields;
my $newcol = defined($ARGV[5]) ? $ARGV[5] : 'SEQUENCE_CONTEXT';
my $header;
my @cols;
my $ffb_cmd;
my $seq;
my $ffb_line;
my ($ref, $alt);

# POSSIBLE_ERROR: Check on anno error, if this is the reason.
$ffb_cmd = "$fastaBinary -fi $refGenome -bed <(vcf2padded_bed.pl $pad $file) -tab -fo stdout";
open(my $ffb_fh, '-|', 'bash', '-c', "$ffb_cmd") || die "Could not open ffb with command $ffb_cmd ($!)";
open(my $snv_fh, "$file") || die "Could not open $file ($!)";

while ($header = <$snv_fh>) {
    last if ($header =~ /^\#CHROM/);
    print $header;
}
die "Error reading $file: $!\n" if $snv_fh->error;  # issue #43
chomp $header;
@cols = split(/\t/, $header);
my %cols = map {$_ => 1} @cols;
push(@cols, $newcol) if (!$cols{$newcol});
say join "\t", @cols;

while (<$snv_fh>) {
    chomp;
    @fields{@cols} = split(/\t/);
    $ffb_line = <$ffb_fh>;
    $ffb_line //= ' ';
    chomp $ffb_line;
    $seq = uc((split(/\t/, $ffb_line))[1]);
    $fields{$newcol} = join ',', substr($seq, 0, $pad), substr($seq, -$pad);
    say join "\t", @fields{@cols};
}
die "Error reading $file: $!\n" if $snv_fh->error;  # issue #43
close $snv_fh;
close($ffb_fh) || die "Could not close fastaFromBed filehandle (broken pipe?) ($!)";