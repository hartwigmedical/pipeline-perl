#!/usr/bin/env perl

# approximate discipline.pm rather than copying to runtime directory
use 5.016_000;
use strictures 2;
no indirect 'fatal';
no multidimensional;
no bareword::filehandles;


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

        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split /\t/, $line;
        my %info_map = map { my ($key, $value) = split /=/; } split /;/, $info;

        # some tools (IGV) will complain if END >= POS.
        # swap around the coordinates of these SVs and
        # (unconventially) report them in reverse.
        if ($info_map{END} >= $pos) {
            say $out_vcf_fh $line;
        } else {
            if ($info_map{CT} eq "5to3") {
                $info_map{CT} = "3to5";
            } elsif ($info_map{CT} eq "3to5") {
                $info_map{CT} = "5to3";
            }

            $info_map{CONSENSUS} = reverse($info_map{CONSENSUS} =~ tr/acgtACGT/tgcaTGCA/r) if exists $info_map{CONSENSUS};
            ($chrom, $info_map{CHR2}) = ($info_map{CHR2}, $chrom);
            ($pos, $info_map{END}) = ($info_map{END}, $pos);

            $info = join ";", map { "$_" . (defined $info_map{$_} ? "=$info_map{$_}" : "") } keys %info_map;
            say $out_vcf_fh join("\t", $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, join("\t", @samples));
        }
    }
    return;
}
