#!/usr/bin/env perl

use discipline;

use File::Spec::Functions;
use File::Temp;
use File::Touch;
use Test::Fatal;
use Test::More;
use Test::Output;
use Test::MockModule;

use lib "t";
# bug in Test::Prereq 1.x needs filename for test dependencies
require "Util.pm"; ## no critic (Modules::RequireBarewordIncludes)

use HMF::Pipeline::Job;


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

sub testDoneFile {
    my ($name, $step, $done_name, $is_reported, $skip) = @_;

    my $temp_dir = File::Temp->newdir();
    my $done_file = catfile($temp_dir, $done_name);
    touch $done_file if $skip;

    my $opt = {};
    my $returned_done_file;
    stdout_is {
        $returned_done_file = HMF::Pipeline::Job::checkDoneFile($name, $step, $is_reported, {log => $temp_dir, mapping => $temp_dir}, $opt)
    }
    $skip ? "WARNING: $done_file exists, skipping\n" : "", $skip ? "notifies of skip" : "no output when not skipping";
    is($returned_done_file, $skip ? undef : $done_file, $skip ? "returns undef when skipping job" : "returns standardised job name when not skipping");
    is_deeply($opt, $is_reported ? {DONE_FILES => [$done_file]} : {}, $is_reported ? "records done file for finalize" : "does not record unreported jobs");

    return;
}

# File::Spec::case_tolerant() is broken
sub case_tolerant {
    my $temp_dir = File::Temp->newdir();
    my $test_file = catfile($temp_dir, "aaa");
    touch $test_file;
    return -f uc $test_file;
}

foreach my $is_reported (0, 1) {
    # mostly backwards-compatibility testing for smooth re-running of parts of older samples
    # tests some normal behaviour at the same time (but not all files - integration tests cover current INIs)
SKIP: {
        skip "case-insensitive filesystem (macOS) always matches some names", 6 if case_tolerant();

        testDoneFile("Freec", undef, "freec.done", $is_reported, 1);
        testDoneFile("QDNAseq", undef, "qdnaseq.done", $is_reported, 1);
        testDoneFile("Strelka", undef, "strelka.done", $is_reported, 1);
    }

    testDoneFile("PerLaneConvert", "core_name", "core_name.done", $is_reported, 1);
    testDoneFile("Map", "core_name", "core_name_bwa.done", $is_reported, 1);
    testDoneFile("PreStats", "core_name", "PreStats_core.done", $is_reported, 1);
    testDoneFile("GermlineCalling", undef, "GermlineCaller.done", $is_reported, 1);
    testDoneFile("GermlineCalling", undef, "VariantCaller.done", $is_reported, 1);
    testDoneFile("GermlineFiltering", undef, "GermlineFilter.done", $is_reported, 1);
    testDoneFile("GermlineFiltering", undef, "VariantFilter.done", $is_reported, 1);
    testDoneFile("GermlineAnnotation", undef, "VariantAnnotation.done", $is_reported, 1);

    testDoneFile("name", "step", "name_step.done", $is_reported, 0);
}

my $temp_dir = File::Temp->newdir();
my $opt = {};
my $returned_done_file;
stdout_is {
    $returned_done_file = HMF::Pipeline::Job::checkReportedDoneFile("name", "step", {log => $temp_dir}, $opt)
}
"", "no output when no done file";
is($returned_done_file, catfile($temp_dir, "name_step.done"), "returns standard done file");
is_deeply($opt, {DONE_FILES => [ catfile($temp_dir, "name_step.done") ]}, "records standard done file for finalize");


is(HMF::Pipeline::Job::fullname("name", undef), "name", "names job with one step");
is(HMF::Pipeline::Job::fullname("name", "step"), "name_step", "names job with multiple steps");


my %ids = map { HMF::Pipeline::Job::getId() => 1 } (1 .. 1000);
is(keys %ids, 1000, "IDs somewhat unique");


my $param = HMF::Pipeline::Job::hold_jid([]);
is($param, "", "no hold_jid param with no jobs");
$param = HMF::Pipeline::Job::hold_jid(["jobid"]);
is($param, "-hold_jid jobid", "hold_jid param with one job");
$param = HMF::Pipeline::Job::hold_jid([ "jobid1", "jobid2" ]);
is($param, "-hold_jid jobid1,jobid2", "hold_jid param with multiple jobs");


stdout_like {
    HMF::Pipeline::Job::submit("qsub", "job_name", "job_id", [], "bash_file.sh", {log => "logs"})
}
qr/Your job [0-9]+ \("job_id"\) has been submitted/, "calls qsub";

my $unwriteable_file = File::Temp->new();
chmod 0000, $unwriteable_file->filename;
my $module = Test::MockModule->new('File::Temp');
$module->mock(new => sub { return $unwriteable_file->filename; });
my $exception = exception { HMF::Pipeline::Job::submit("qsub", "job_name", "job_id", [], "bash_file.sh", {log => "logs"}) };
like($exception, qr/failed to write qsub parameter file $unwriteable_file/, "fails when qsub options cannot be written");


sub testJobFromTemplate {
    my ($name, $step, $is_reported_job, $hold_jids, $expected_job_name, $expected_hold_jids, $skip) = @_;

    my $job_dir = File::Temp->newdir();
    my $log_dir = File::Temp->newdir();

    my (@check_args, @template_args, @submit_args);
    my $module = Test::MockModule->new('HMF::Pipeline::Job');
    my $expected_job_id = "${expected_job_name}_randomID";
    my $expected_done_file = catfile($log_dir, "${expected_job_name}.done");
    $module->mock(getId => sub { return "randomID" });
    $module->mock(checkDoneFile => sub { @check_args = @_; return $skip ? undef : $expected_done_file });
    $module->mock(writeFromTemplate => sub { @template_args = @_ });
    $module->mock(submit => sub { @submit_args = @_ });

    my $qsub = "qsub param1 param2";
    my $dirs = {job => $job_dir, log => $log_dir};
    my $opt = {};
    my $job_id = HMF::Pipeline::Job::fromTemplate($name, $step, $is_reported_job, $qsub, $hold_jids, $dirs, $opt);
    is_deeply([@check_args], [ $name, $step, $is_reported_job, $dirs, $opt ], "checked .done file",);
    is_deeply(
        [@template_args],
        $skip
        ? []
        : [
            "${name}.sh.tt",
            catfile($job_dir, "${expected_job_id}.sh"),
            done_file => $expected_done_file,
            dirs => $dirs,
            opt => $opt,
        ],
        $skip ? "template write skipped for skipped job" : "wrote template"
    );
    #<<< no perltidy
    is_deeply(
        [@submit_args],
        $skip
        ? []
        : [
            $qsub,
            $expected_job_name,
            $expected_job_id,
            $expected_hold_jids,
            catfile($job_dir, "${expected_job_id}.sh"),
            $dirs,
        ],
        $skip ? "submission skipped for skipped job" : "submitted job"
    );
    #>>> no perltidy
    is($job_id, $skip ? undef : $expected_job_id, $skip ? "no job id for skipped job" : "job id returned");

    return;
}

foreach my $skip (0, 1) {
    foreach my $is_reported (0, 1) {
        testJobFromTemplate("name", "step", $is_reported, [], "name_step", [], $skip);
        testJobFromTemplate("name", undef, $is_reported, [], "name", [], $skip);

        testJobFromTemplate("name", "step", $is_reported, ["jobid"], "name_step", ["jobid"], $skip);
        testJobFromTemplate("name", undef, $is_reported, ["jobid"], "name", ["jobid"], $skip);
        testJobFromTemplate("name", undef, $is_reported, [undef], "name", [], $skip);
        testJobFromTemplate("name", undef, $is_reported, [ "jobid1", undef, "jobid2" ], "name", [ "jobid1", "jobid2" ], $skip);
    }
}


my (@job_args, @qsub_args);
$module = Test::MockModule->new('HMF::Pipeline::Job');
$module->mock(fromTemplate => sub { @job_args = @_; return "jobid"; });
$module->mock(qsubSimple => sub { @qsub_args = @_ });

is(HMF::Pipeline::Job::markDone(undef, [], {}, {}), undef, "skips marking done file already detected as written");
is_deeply([@job_args], [], "no job run");

my $job_id = HMF::Pipeline::Job::markDone(catfile($temp_dir, "something.done"), ["jobid"], {}, {});
is($job_id, "jobid", "runs mark done job");

done_testing();
