#!/usr/bin/env perl

use strict;
use warnings;

use File::Spec::Functions;
use File::Temp;
use Test::Dir;
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


$path = catfile("some", "directory", "sample_flowcell_index_lane_R1_suffix.fastq.gz");
$opt = {FASTQ => {$path => 1}};
addSamples($opt);
is_deeply(
    $opt, {
        FASTQ => {$path => 1},
        SAMPLES => {sample => [$path]},
        RUNNING_JOBS => {sample => []},
    },
    "sets up FASTQ sample"
);

my $other_path = catfile("some", "directory", "sample_flowcell_index_lane_R2_suffix.fastq.gz");
$opt = {FASTQ => {$path => 1, $other_path => 1}};
addSamples($opt);
is_deeply(
    $opt, {
        FASTQ => {$path => 1, $other_path => 1},
        SAMPLES => {sample => [ $path, $other_path ]},
        RUNNING_JOBS => {sample => []},
    },
    "sets up multiple FASTQ samples"
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
            SAMPLES => {empty => [$path]},
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

my $temp_dir = File::Temp->newdir();
my ($test_path_a, $test_path_b);

$test_path_a = catfile($temp_dir, "dir_a");
HMF::Pipeline::Config::makePaths($test_path_a);
dir_exists_ok($test_path_a, "makes single path");

HMF::Pipeline::Config::makePaths($test_path_a);
dir_exists_ok($test_path_a, "makes pre-existing path");

$test_path_a = catfile($temp_dir, "dir_a");
$test_path_a = catfile($temp_dir, "dir_b");
HMF::Pipeline::Config::makePaths($test_path_a, $test_path_b);
dir_exists_ok($test_path_a, "makes first path");
dir_exists_ok($test_path_a, "makes second path");

$test_path_a = catfile($temp_dir, "nested", "dir_a");
HMF::Pipeline::Config::makePaths($test_path_a);
dir_exists_ok($test_path_a, "makes nested path");

$temp_dir = File::Temp->newdir();
chmod 0000, $temp_dir;
$test_path_a = catfile($temp_dir, "dir_a");
$exception = exception { HMF::Pipeline::Config::makePaths($test_path_a) };
dir_not_exists_ok($test_path_a, "no directory on failure");
like($exception, qr(Couldn't create directories: .*/dir_a: Permission denied), "fails when directory cannot be created");

$temp_dir = File::Temp->newdir();
my $output_dir = catfile($temp_dir, "module");

my $dirs = HMF::Pipeline::Config::createDirs($output_dir, extra => "extra_dir", nested => catfile("base", "nested_dir"));
is_deeply(
    $dirs, {
        out => $output_dir,
        tmp => catfile($output_dir, "tmp"),
        log => catfile($output_dir, "logs"),
        job => catfile($output_dir, "jobs"),
        extra => catfile($output_dir, "extra_dir"),
        nested => catfile($output_dir, "base", "nested_dir"),
    },
    "returns standard and extra directory mapping"
);
dir_exists_ok($output_dir, "makes out dir");
dir_exists_ok(catfile($output_dir, "tmp"), "makes tmp dir");
dir_exists_ok(catfile($output_dir, "logs"), "makes log dir");
dir_exists_ok(catfile($output_dir, "jobs"), "makes job dir");
dir_exists_ok(catfile($output_dir, "extra_dir"), "makes extra dir");
dir_exists_ok(catfile($output_dir, "base", "nested_dir"), "makes nested extra dir");

my $subdir;

$subdir = HMF::Pipeline::Config::addSubDir($dirs, "subdir");
is($subdir, catfile($output_dir, "subdir"), "returns sub-directory");
dir_exists_ok(catfile($output_dir, "subdir"), "adds sub-directory");

$subdir = HMF::Pipeline::Config::addSubDir($dirs, "subdir");
is($subdir, catfile($output_dir, "subdir"), "returns already-existing sub-directory");
dir_exists_ok(catfile($output_dir, "subdir"), "adds already-existing sub-directory");

$subdir = HMF::Pipeline::Config::addSubDir($dirs, catfile("base", "subdir"));
is($subdir, catfile($output_dir, "base", "subdir"), "returns nested sub-directory");
dir_exists_ok(catfile($output_dir, "base", "subdir"), "adds nested sub-directory");

my $chrs = HMF::Pipeline::Config::getChromosomes({GENOME => catfile("t", "data", "empty.fasta")});
is_deeply($chrs, [1], "gets chromosomes");

done_testing();
