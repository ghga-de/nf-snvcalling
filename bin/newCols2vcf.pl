#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
### This script adds columns to a vcf file (or updates them if columns with these names are already present).
### Both files must be sorted in the same order (including chromosomal order).
### It is for exact matches, e.g. to merge annovar output to a vcf.
### It is crucial that there is no entry in the newcol File which does not match an entry in the vcf file
### Missing entries in the newcol file are fine
### Chromosomal identifiers in the newcol file may be plain numbers while longer identifiers are used in the vcf file e.g. after annovar; use chrPrefix and chrSuffix then
### (long chr identifiers in the newcol file and short identifiers in the vcf file are not supported)
# editted by kuebra.narci@dkfz.de to correct error handling - 23.10.2024 

use strict;
use warnings;
use v5.10;


my %opts;
BEGIN {
  use Getopt::Long;
    

  %opts = ('aFileType' => 'vcf',
	   'chrPrefix' => '',
	   'chrSuffix' => '',
	   'endAdd' => 0,
	   'aChrPosEnd' => '',
	   'bChrPosEnd' => '',
	  );
  
  GetOptions( 'vcfFile=s' => \$opts{vcf},
              'newColFile=s' => \$opts{newcolfile},
	      'newColHeader=s' => \$opts{newcolheader},
	      'chrPrefix:s' => \$opts{chrPrefix},
	      'chrSuffix:s' => \$opts{chrSuffix},
	      'aChrPosEnd:s' => \$opts{aChrPosEnd},
	      'bChrPosEnd=s' => \$opts{bChrPosEnd},
              'endAdd:i' => \$opts{endAdd},
              'reportColumns=s' => \$opts{reportColumns},
	      'aFileType:s' => \$opts{aFileType},
            );

  


  my @bCPE = split(',',$opts{bChrPosEnd});
  if (@bCPE < 3) {
    die "Error: bChrPosEnd not specified correctly";
  }


  require constant;
  constant->import(CHRPREFIX => $opts{chrPrefix});
  constant->import(CHRSUFFIX => $opts{chrSuffix});
  constant->import(AFILETYPE => $opts{aFileType});
  constant->import(ENDADD => $opts{endAdd});
  constant->import(BCHRCOL => $bCPE[0]);
  constant->import(BPOSCOL => $bCPE[1]);
  constant->import(BENDCOL => $bCPE[2]);
  
  my @aCPE;
  # Check if the required options are defined
  if (!defined $opts{aFileType}) {
      die "Error: 'aFileType' is not specified.\n";
  }

  if ($opts{aFileType} eq 'custom') {
      if (!defined $opts{aChrPosEnd}) {
          die "Error: 'aChrPosEnd' must be specified for custom file types.\n";
      }
      @aCPE = split(',', $opts{aChrPosEnd});
      die "Invalid specification of aChrPosEnd: at least 3 elements are required." if (@aCPE < 3);
  } elsif ($opts{aFileType} eq 'vcf') {
      @aCPE = (0, 1, 'NA');
  } else {
      die "Invalid a-file type specified: '$opts{aFileType}'. Please use 'custom' or 'vcf'.\n";
  }


  constant->import(ACHRCOL => $aCPE[0]);
  constant->import(APOSCOL => $aCPE[1]);
  constant->import(AENDCOL => $aCPE[2]);
}

my @repCols = split(',',$opts{reportColumns});
(@repCols) || die "No report columns set";
my @newcols = split(',',$opts{newcolheader});


#  @nc_name{@repCols} = @newcols;

open(my $vcf_fh, "<", "$opts{vcf}") || die "Could not open VCF file $opts{vcf}";
open(my $nc_fh, "<", "$opts{newcolfile}") || die "Could not open new column file $opts{newcolfile}";

my $header;
while ($header = <$vcf_fh>) {
    last if ($header =~ /^\#CHR/); # Identify the header line
    print $header; # Print metadata
}
chomp($header);
my @columns = split(/\t/, $header);
my @ori_columns = @columns;

foreach my $new_col (@newcols) {
    unless ($header =~ /\Q$new_col\E/) {
        push(@columns, $new_col);
        $header .= "\t$new_col";
    }
}
say $header;
my %f1_hash;

my ($l1, $l2, $end);
my @f2_fields;
NC_LOOP: while ($l2=<$nc_fh>) {
  next if ($l2 =~ /^\#/);
  chomp $l2;
  @f2_fields = split(/\t/, $l2);
  while ($l1=<$vcf_fh>) {
    chomp $l1;
    @f1_hash{@ori_columns} = split(/\t/, $l1);
    if (AFILETYPE() eq 'vcf') {
      if ($f1_hash{INFO} =~ /END=(\d+)/) {
	$end = $1;
      } else {
	$end = $f1_hash{POS} + length($f1_hash{REF}) - 1;
      }
    } else {
      $end = $f1_hash{$columns[AENDCOL()]};
    }
    if (CHRPREFIX().$f2_fields[BCHRCOL()].CHRSUFFIX() eq $f1_hash{@columns[ACHRCOL()]} && $f2_fields[BPOSCOL()] == $f1_hash{@columns[APOSCOL()]} && $f2_fields[BENDCOL()]+ENDADD() == $end) {
      @f1_hash{@newcols} = @f2_fields[@repCols];
      say join "\t", @f1_hash{@columns};
      next NC_LOOP;
    } else {
      @f1_hash{@newcols} = ('.') x @newcols;
      say join "\t", @f1_hash{@columns};
      # warn "No new column line found for vcf line $l1";
    }
  }
  warn "Line $l2 left over in new column file";
}
#when I am here the NC file has ended; write out every remaining line from VCF
while ($l1=<$vcf_fh>) {
  chomp $l1;
  @f1_hash{@columns} = split(/\t/, $l1);
  @f1_hash{@newcols} = ('.') x @newcols;
  say join "\t", @f1_hash{@columns};
}
close $vcf_fh or die "Could not close VCF file $opts{vcf}";
close $nc_fh or die "Could not close new col file $opts{newcolfile}";