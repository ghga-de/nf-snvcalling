#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

# Fixed 2026-08-12 @kubranarci: I/O error detection (issue #43): autodie for open/close, and an explicit
# read-error check after each read loop (<FH> returns undef on a read error just like at EOF)

use strict;
use warnings;
use autodie;
use IO::Handle;

my @bins;
my $max = 50000;
my $median;
my $outfile = $ARGV[1];

my $in_fh;
if (!defined $ARGV[0] || $ARGV[0] eq '-') {
    $in_fh = *STDIN;
} else {
    open $in_fh, '<', $ARGV[0];
}

my $count=0;
my @head;
while (<$in_fh>) {
	print $_;
	if($_ =~ /^#CHR/){
		@head = split("\t", $_);
		last;
	}
}
die "Error reading input: $!\n" if $in_fh->error;  # issue #43

my $j = 0;
my $dpcol;
foreach(@head){
	if($_ =~ /^INFO_control/){
		$dpcol = $j;
	}
	$j++;
}

while (<$in_fh>) {
	print $_;
	my @l = split("\t", $_);
	my ($dp) = $l[$dpcol] =~ /^DP=(\d+);/;
	next if($dp == 0);
	if($dp <= $max){
		$bins[$dp]++;
	}else{
		$bins[$max]++;
	}
	$count++;
}
die "Error reading input: $!\n" if $in_fh->error;  # issue #43

my $mid = $count/2;

my $medcount = 0;
my $i = 0;
while($i < @bins){
	if(!defined $bins[$i]){
		$bins[$i] = 0;
	}
	$medcount += $bins[$i];
	last if($medcount >= $mid);
	$i++;
}

$median = 5*$i;
close $in_fh;
open my $out_fh, '>', $outfile;
print $out_fh $median, "\n";
close $out_fh;
