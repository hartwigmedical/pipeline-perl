#!/usr/bin/env perl

use strict;
use warnings;

use File::Spec::Functions;
use File::Temp;
use File::Touch;
use Test::More;
use Test::Output;

use lib "t";
# bug in Test::Prereq 1.x needs filename for test dependencies
require "Util.pm"; ## no critic (Modules::RequireBarewordIncludes)

use HMF::Pipeline::Job;


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

sub testDoneFileSkipped {
    my ($name, $step, $done_name) = @_;

    my $temp_dir = File::Temp->newdir();
    my $done_file = catfile($temp_dir, $done_name);
    touch $done_file;

    my $opt = {};
    my $returned_done_file;
    stdout_is {
        $returned_done_file = HMF::Pipeline::Job::checkDoneFile($name, $step, 1, {log => $temp_dir, mapping => $temp_dir}, $opt)
    }
    "WARNING: $done_file exists, skipping\n", "notifies of skip";
    is($returned_done_file, undef, "returns undef when skipping job");
    is_deeply($opt, {DONE_FILES => [$done_file]}, "records done file for finalize");

    return;
}

# File::Spec::case_tolerant() is broken
sub case_tolerant {
    my $temp_dir = File::Temp->newdir();
    my $test_file = catfile($temp_dir, "aaa");
    touch $test_file;
    return -f uc $test_file;
}

# backwards-compatibility testing for smooth re-running of parts of older samples
# covers normal skipping at the same time
SKIP: {
    skip "case-insensitive filesystem (macOS) always matches some names", 6 if case_tolerant();

    testDoneFileSkipped("Freec", undef, "freec.done");
    testDoneFileSkipped("QDNAseq", undef, "qdnaseq.done");
    testDoneFileSkipped("Strelka", undef, "strelka.done");
    testDoneFileSkipped("Varscan", undef, "varscan.done");
    testDoneFileSkipped("Freebayes", undef, "freebayes.done");
    testDoneFileSkipped("Mutect", undef, "mutect.done");
}

testDoneFileSkipped("PerLaneConvert", "core_name", "core_name.done");
testDoneFileSkipped("Map", "core_name", "core_name_bwa.done");
testDoneFileSkipped("PreStats", "core_name", "PreStats_core.done");
testDoneFileSkipped("GermlineCalling", undef, "GermlineCaller.done");
testDoneFileSkipped("GermlineCalling", undef, "VariantCaller.done");
testDoneFileSkipped("GermlineFiltering", undef, "GermlineFilter.done");
testDoneFileSkipped("GermlineFiltering", undef, "VariantFilter.done");
testDoneFileSkipped("GermlineAnnotation", undef, "VariantAnnotation.done");


my $temp_dir = File::Temp->newdir();
my $opt = {};
my $returned_done_file;
stdout_is {
    $returned_done_file = HMF::Pipeline::Job::checkDoneFile("name", "step", 1, {log => $temp_dir}, $opt)
}
"", "no output when no done file";
is($returned_done_file, catfile($temp_dir, "name_step.done"), "returns standard done file");
is_deeply($opt, {DONE_FILES => [ catfile($temp_dir, "name_step.done") ]}, "records standard done file for finalize");

$opt = {};
stdout_is {
    $returned_done_file = HMF::Pipeline::Job::checkReportedDoneFile("name", "step", {log => $temp_dir}, $opt)
}
"", "no output when no done file";
is($returned_done_file, catfile($temp_dir, "name_step.done"), "returns standard done file");
is_deeply($opt, {DONE_FILES => [ catfile($temp_dir, "name_step.done") ]}, "records standard done file for finalize");

$opt = {};
$returned_done_file = HMF::Pipeline::Job::checkDoneFile("name", "step", 0, {log => $temp_dir}, $opt);
is($returned_done_file, catfile($temp_dir, "name_step.done"), "returns standard done file for non-reported job");
is_deeply($opt, {}, "does not record standard done file for non-reported job");


done_testing();
