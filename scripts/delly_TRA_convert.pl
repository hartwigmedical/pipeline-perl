#! /usr/bin/perl

use strict;

my $vcf = shift;
my $vcf_out = $vcf;
$vcf_out =~ s/.vcf$/_CONVERT.vcf/;
open OUT, ">$vcf_out" or die "Cannot create file\n$!\n";

open VCF, $vcf;
while(my $line = <VCF>) {
	chomp($line);
	print OUT $line . "\n" if $line =~ /^#/;
	next if $line =~ /^#/;
	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split/\t/, $line;
	my $chrom2 = $1 if $info =~ /CHR2=(\w+);/;
	my $end = $1 if $info =~ /END=(\d+);/;
	my $ct = $1 if $info =~ /CT=(\w+);/;
	my $consensus = $1 if $info =~ /CONSENSUS=(\w+)/;
	print OUT $line . "\n" if $end >= $pos;
	if ($end < $pos) {
		if ($ct eq "5to3") {
			$ct = "3to5";
		} elsif ($ct eq "3to5") {
			$ct = "5to3";
		} else {
		
		}
		$consensus =~ tr/acgtACGT/tgcaTGCA/;
		$consensus = reverse($consensus);
		$info =~ s/CHR2=(\w+)/CHR2=$chrom/;
		$info =~ s/END=(\d+)/END=$pos/;
		$info =~ s/CT=(\w+)/CT=$ct/;
		$info =~ s/CONSENSUS=(\w+)/CONSENSUS=$consensus/;
		print OUT join("\t", $chrom2, $end, $id, $ref, $alt, $qual, $filter, $info, $format, join("\t", @samples)) . "\n";
	}
}
print "Done.\n";
close VCF;
close OUT;