#!/usr/bin/perl -w
use strict;
use Number::Format;

#collect flagstats and tally up the totals.

my $fn=new Number::Format(	-thousands_sep 	=> '.',
				-decimal_point	=> ',');

my @files = `find \$PWD -iname \"*MERGED*.flagstat\"`;
print scalar @files, " flagstats found\n";

my %totals;

foreach my $file (@files) {
    chomp($file);
    print "working on $file...\n";
    open IN, $file;
    my ($mapped, $dups) = (0,0);
    while (my $line=<IN>) {
	print "\t$line";
	if ($line =~ /(\d+) \+ \d in total/) {
	    $totals{'total'} += $1;
	}elsif($line =~ /(\d+) \+ \d paired in/) {
	    $totals{'paired'} += $1;
	}elsif($line =~ /(\d+) \+ \d properly paired/) {
	    $totals{'ppaired'} += $1;
	}elsif($line =~ /(\d+) \+ \d duplicates/) {
	    $totals{'dups'} += $1;
	    $dups = $1;
	}elsif($line =~ /(\d+) \+ \d mapped \(/) {
	    $totals{'mapped'} += $1;
	    $mapped = $1;
	}
    }
    print "\n\t";
    print (100*$dups/$mapped);
    print " %duplication\n\n";
    
}


print "Total raw reads: ".$fn->format_number($totals{'total'})," reads (", $fn->format_number(100*$totals{'total'})," bp)\n";
print "Total mapped reads: ".$fn->format_number($totals{'mapped'})," reads (", $fn->format_number(100*$totals{'mapped'})," bp)\n";
print "Average mapped per lib: ".$fn->format_number(int($totals{'total'}/(scalar @files)))," reads\n";
print "Average dups per lib: ".$fn->format_number(int($totals{'dups'}/(scalar @files)))," reads\n";
print "Average dups % per lib: ",$fn->format_number( 100*($totals{'dups'}/$totals{'mapped'}))," %\n";