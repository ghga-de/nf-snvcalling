#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/COWorkflowsBasePlugin/LICENSE).
#

use strict;
use warnings;
use v5.10;

my $file = $ARGV[0];
my $refGenome = $ARGV[1];
my $pad = $ARGV[2];
my %fields;
my $newcol = defined($ARGV[3]) ? $ARGV[3] : 'SEQUENCE_CONTEXT';
my $header;
my @cols;
my $bash = '/bin/bash';
my $ffb_cmd;
my $seq;
my $ffb_line;
my ($ref, $alt);

# POSSIBLE_ERROR: Check on anno error, if this is the reason.
$ffb_cmd = "fastaFromBed -fi $refGenome -bed <(perl vcf2padded_bed.pl $pad $file) -tab -fo stdout";
open(my $ffb_fh, '-|', $bash, '-c', "$ffb_cmd") || die "Could not open ffb with command $ffb_cmd ($!)";
open(my $snv_fh, "$file") || die "Could not open $file ($!)";

while ($header = <$snv_fh>) {
    last if ($header =~ /^\#CHROM/);
    print $header;
}
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
close $snv_fh;
close($ffb_fh) || die "Could not close fastaFromBed filehandle (broken pipe?) ($!)";
