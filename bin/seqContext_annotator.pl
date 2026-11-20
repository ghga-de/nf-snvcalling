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

# fastaFromBed skips intervals it cannot extract (contig missing from the FASTA, or
# beyond the contig end), so its output lines cannot be paired with the VCF lines by
# position alone: one skipped interval shifted the contexts of all following variants
# and left the last ones empty (#45). Each fastaFromBed line is therefore matched to
# its variant by its interval name (chr:start-end). Variants without a sequence get
# an all-N context, and contexts at contig borders are padded with N, so every
# context has $pad bases on both sides.
my $fai = "$refGenome.fai";
my %contig_length;
if (-e $fai) {
    open my $fai_fh, '<', $fai;
    while (<$fai_fh>) {
        my ($name, $length) = split(/\t/);
        $contig_length{$name} = $length;
    }
    die "Error reading $fai: $!\n" if $fai_fh->error;
    close $fai_fh;
} else {
    $fai = '';
}

$ffb_cmd = "$fastaBinary -fi $refGenome -bed <(vcf2padded_bed.pl $pad $file $fai) -tab -fo stdout";
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

my $pending;          # fastaFromBed line read ahead, not yet matched to a variant
my $no_sequence = 0;  # variants without a sequence from fastaFromBed
while (<$snv_fh>) {
    chomp;
    @fields{@cols} = split(/\t/);
    my ($chr, $pos, $ref_allele) = @fields{'#CHROM', 'POS', 'REF'};

    # interval as written by vcf2padded_bed.pl (0-based start, end clipped to the contig)
    my $start = ($pos - $pad - 1 > 0) ? $pos - $pad - 1 : 0;
    my $end   = $pos + length($ref_allele) + $pad - 1;
    $end = $contig_length{$chr} if (exists $contig_length{$chr} && $end > $contig_length{$chr});
    my $interval = "$chr:$start-$end";

    if (!defined $pending) {
        $pending = <$ffb_fh>;
        chomp $pending if (defined $pending);
    }
    my ($before, $after) = ('', '');
    if (defined $pending && (split(/\t/, $pending))[0] eq $interval) {
        $seq = uc((split(/\t/, $pending))[1]);
        $before = substr($seq, 0, $pos - 1 - $start);
        $after  = substr($seq, $pos - 1 - $start + length($ref_allele));
        undef $pending;
    } else {
        $no_sequence++;
    }
    # always $pad bases on each side; N where no sequence exists
    $before = substr(('N' x $pad) . $before, -$pad);
    $after  = substr($after . ('N' x $pad), 0, $pad);
    $fields{$newcol} = join ',', $before, $after;
    say join "\t", @fields{@cols};
}
die "Error reading $file: $!\n" if $snv_fh->error;  # issue #43
die "fastaFromBed returned sequences that match no variant (first: $pending)\n" if (defined $pending);
my $leftover = <$ffb_fh>;
die "fastaFromBed returned sequences that match no variant (first: $leftover)" if (defined $leftover);
warn "No sequence context for $no_sequence variant(s) (contig not in the reference?); their $newcol is all N\n" if ($no_sequence);
close $snv_fh;
close($ffb_fh) || die "Could not close fastaFromBed filehandle (broken pipe?) ($!)";
