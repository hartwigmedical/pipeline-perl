#!/usr/bin/env perl

use strict;
use warnings;

use File::Spec::Functions;
use Test::Fatal;
use Test::More;

use lib "t";
# bug in Test::Prereq 1.x needs filename for test dependencies
require "Util.pm"; ## no critic (Modules::RequireBarewordIncludes)

use HMF::Pipeline::Config qw(addSamples);


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

sub testParse {
    my ($string) = @_;

    my $opt = {};
    open my $fh, "<", $string;
    HMF::Pipeline::Config::parseFile($fh, $opt);
    close $fh;
    return $opt;
}

my $opt;
my $exception;
my $path;

$opt = {};
$opt = testParse(\<<"_EOF_");
VALID	value
_EOF_
is_deeply($opt, {VALID => "value"}, "parses simple key and value");

$opt = {};
$exception = exception {
    $opt = testParse(\<<"_EOF_"); };
INVALID
_EOF_
like($exception, qr/\QKey 'INVALID' is missing a value - is it badly formatted?\E/, "refuses key with no value");

$opt = {};
$exception = exception {
    $opt = testParse(\<<"_EOF_") };
INVALID value
_EOF_
like($exception, qr/\QKey 'INVALID value' is missing a value - is it badly formatted?\E/, "refuses key with space separator");

$opt = {};
$opt = testParse(\<<"_EOF_");
FASTQ	file1
FASTQ	file2
_EOF_
is_deeply($opt, {FASTQ => {file1 => 1, file2 => 1}}, "parses FASTQs into set");

$opt = {};
$opt = testParse(\<<"_EOF_");
BAM	file1
BAM	file2
_EOF_
is_deeply($opt, {BAM => {file1 => 1, file2 => 1}}, "parses BAMs into set");


my $test_ini_path = catfile("t", "data", "test.ini");
$path = catfile(HMF::Pipeline::Config::pipelinePath(), $test_ini_path);
$opt = {};
$opt = testParse(\<<"_EOF_");
INIFILE	$test_ini_path
_EOF_
is_deeply($opt, {INIFILE => [$path], VALID => "other_value"}, "parses INIFILE");

$path = catfile(HMF::Pipeline::Config::pipelinePath(), $test_ini_path);
$opt = {};
$opt = testParse(\<<"_EOF_");
VALID	value
INIFILE	$test_ini_path
_EOF_
is_deeply($opt, {INIFILE => [$path], VALID => "other_value"}, "later INIFILE overrides earlier key");

$path = catfile(HMF::Pipeline::Config::pipelinePath(), $test_ini_path);
$opt = {};
$opt = testParse(\<<"_EOF_");
INIFILE	$test_ini_path
VALID	value
_EOF_
is_deeply($opt, {INIFILE => [$path], VALID => "value"}, "later key overrides earlier INIFILE");


$path = catfile("some", "directory", "sample_something_something.fastq.gz");
$opt = {FASTQ => {$path => 1}};
addSamples($opt);
is_deeply(
    $opt, {
        FASTQ => {$path => 1},
        SAMPLES => {sample => $path},
        RUNNING_JOBS => {sample => []},
    },
    "sets up FASTQ sample"
);

SKIP: {
    skip "no SAMTOOLS_PATH set", 2 if not $ENV{SAMTOOLS_PATH};

    $path = Util::convertToBam(catfile("t", "data", "empty.sam"));
    $opt = {
        BAM => {$path => 1},
        SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
        GENOME => catfile("t", "data", "empty.fasta"),
    };
    addSamples($opt);
    is_deeply(
        $opt, {
            BAM => {$path => 1},
            SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
            GENOME => catfile("t", "data", "empty.fasta"),
            SAMPLES => {empty => $path},
            RUNNING_JOBS => {empty => []},
        },
        "sets up BAM sample"
    ) or diag explain $opt;


    my $other_path = Util::convertToBam(catfile("t", "data", "other_empty.sam"));
    $opt = {
        BAM => {$path => 1, $other_path => 1},
        SAMTOOLS_PATH => $ENV{SAMTOOLS_PATH},
        GENOME => catfile("t", "data", "empty.fasta"),
    };
    $exception = exception { HMF::Pipeline::Config::addSamples($opt) };
    like($exception, qr/sample 'empty' from $other_path already used by $path/, "refuses duplicate BAM sample") or diag explain $opt;
}

done_testing();
