#!/usr/bin/env perl

# TODO: use discipline.pm
use 5.016_000;
use strict;
use warnings;


say "Starting TRA conversion";
my $in_vcf = shift;
my $out_vcf = shift;
open my $out_vcf_fh, ">", $out_vcf or die "Cannot create file $out_vcf: $!";
open my $in_vcf_fh, "<", $in_vcf or die "Cannot read file $in_vcf: $!";
convertTra($in_vcf_fh, $out_vcf_fh);
say "Done.";
close $in_vcf_fh;
close $out_vcf_fh;


sub convertTra {
    while (my $line = <$in_vcf_fh>) {
        chomp($line);
        say $out_vcf_fh $line and next if $line =~ /^#/;

        # TODO: consider parsing properly e.g. using PyVCF
        # this has potential bugs with INFO field ordering, substring matches etc.
        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split /\t/, $line;
        (my $chrom2 = $info) =~ s/CHR2=(\w+);/$1/;
        (my $end = $info) =~ s/END=(\d+);/$1/;
        (my $ct = $info) =~ s/CT=(\w+);/$1/;
        (my $consensus = $info) =~ s/CONSENSUS=(\w+)/$1/;

        # TODO: bug? need chromosomes to match also?
        if ($end >= $pos) {
            say $out_vcf_fh $line;
        } else {
            if ($ct eq "5to3") {
                $ct = "3to5";
            } elsif ($ct eq "3to5") {
                $ct = "5to3";
            } else {
                # TODO: what about 3to3 and 5to5?
                die "ERROR: unrecognised connection type CT=$ct";
            }
            $consensus =~ tr/acgtACGT/tgcaTGCA/;
            $consensus = reverse($consensus);
            $info =~ s/CHR2=(\w+)/CHR2=$chrom/;
            $info =~ s/END=(\d+)/END=$pos/;
            $info =~ s/CT=(\w+)/CT=$ct/;
            $info =~ s/CONSENSUS=(\w+)/CONSENSUS=$consensus/;
            say $out_vcf_fh join("\t", $chrom2, $end, $id, $ref, $alt, $qual, $filter, $info, $format, join("\t", @samples));
        }
    }
    return;
}
