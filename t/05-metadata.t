#!/usr/bin/env perl

use discipline;

use File::Basename;
use File::Path qw(make_path);
use File::Spec::Functions;
use File::Temp;
use Test::Fatal;
use Test::More;

use HMF::Pipeline::Functions::Metadata qw(parse metaSampleName sampleControlNames linkArtefact linkVcfArtefacts linkBamArtefacts readLinks writeLinks);


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

my $temp_dir = File::Temp->newdir();
my $opt = {OUTPUT_DIR => $temp_dir, RUN_NAME => basename($temp_dir)};

my $metadata_path = catfile($temp_dir, "metadata");
HMF::Pipeline::Functions::Metadata::writeJson(
    $metadata_path, {
        type_a => "sample_a",
        type_b => "sample_b",
        random_unused_key => undef,
    }
);

my $metadata = parse($opt);
is_deeply(
    $metadata, {
        type_a => "sample_a",
        type_b => "sample_b",
        random_unused_key => undef,
    },
    "reads metadata"
);

is(metaSampleName("sample_a", $opt), "type_a", "names sample file from metadata");
is(metaSampleName("sample_b", $opt), "type_b", "names another sample file from metadata");
is(metaSampleName("sample_c", $opt), "sample", "default name for missing sample file");

$opt->{SAMPLES} = {"sample_a" => "orig_a.bam", "sample_b" => "orig_b.bam"};
$opt->{BAM_FILES} = {"sample_a" => "filename_a.bam", "sample_b" => "filename_b.bam"};
$opt->{RUNNING_JOBS} = {};
linkBamArtefacts($opt);
ok(exists $opt->{LINKS}, "links stored in \$opt");
is_deeply(
    $opt->{LINKS}, {
        type_a_bam => catfile("sample_a", "mapping", "filename_a.bam"),
        type_a_bai => catfile("sample_a", "mapping", "filename_a.bam.bai"),
        type_b_bam => catfile("sample_b", "mapping", "filename_b.bam"),
        type_b_bai => catfile("sample_b", "mapping", "filename_b.bam.bai"),
    },
    "bam filename links stored"
);

delete $opt->{LINKS};
linkArtefact("filename_a", "artefact_a", $opt);
linkArtefact("filename_b", "artefact_b", $opt);
ok(exists $opt->{LINKS}, "links stored in \$opt");
is_deeply($opt->{LINKS}, {artefact_a => "filename_a", artefact_b => "filename_b"}, "artefact filename links stored");

delete $opt->{LINKS};
linkVcfArtefacts("filename", "artefact", $opt);
ok(exists $opt->{LINKS}, "vcf links stored in \$opt");
is_deeply($opt->{LINKS}, {artefact_vcf => "filename", artefact_vcf_index => "filename.tbi"}, "vcf artefact links stored");

my $links = readLinks($opt);
is_deeply($links, {}, "empty links when nothing written to file");

make_path(catfile($temp_dir, "logs"));
delete $opt->{LINKS};
linkArtefact("filename_a", "artefact_a", $opt);
linkArtefact("filename_b", "artefact_b", $opt);
writeLinks($opt);
my $links_path = catfile($temp_dir, "logs", "links.json");
$links = HMF::Pipeline::Functions::Metadata::readJson($links_path);
is_deeply(
    $links, {
        artefact_a => "filename_a",
        artefact_b => "filename_b",
    },
    "links written to file"
);

is_deeply(readLinks($opt), $links, "reads same links when written to file");

$opt->{LINKS} = {"artefact_a" => catfile($temp_dir, "filename_a")};
writeLinks($opt);
is_deeply(readLinks($opt), {"artefact_a" => "filename_a"}, "converts absolute path in existing links.json");

HMF::Pipeline::Functions::Metadata::writeJson(
    $metadata_path, {
        ref_sample => "ref",
        tumor_sample => "tumor",
    }
);
my $names = [ sampleControlNames($opt) ];
is_deeply($names, [ "ref", "tumor", "ref_tumor" ], "gets sample/control/joint names");

my $exception;
HMF::Pipeline::Functions::Metadata::writeJson(
    $metadata_path, {
        tumor_sample => "tumor",
    }
);
$exception = exception { my $names = [ sampleControlNames($opt) ] };
like($exception, qr/metadata missing ref_sample/, "detects missing ref sample");

HMF::Pipeline::Functions::Metadata::writeJson(
    $metadata_path, {
        ref_sample => "ref",
    }
);
$exception = exception { my $names = [ sampleControlNames($opt) ] };
like($exception, qr/metadata missing tumor_sample/, "detects missing tumor sample");

$temp_dir->DESTROY();
$exception = exception { $metadata = HMF::Pipeline::Functions::Metadata::readJson($links_path) };
like($exception, qr/Can't open $links_path:/, "fails to open missing file");
$exception = exception { HMF::Pipeline::Functions::Metadata::writeJson($links_path, {}) };
like($exception, qr/Can't open $links_path:/, "fails to write to non-existent directory");

$exception = exception { $metadata = parse($opt) };
like($exception, qr/Can't open $metadata_path:/, "error when metadata missing");

$exception = exception { $metadata = parse($opt, {required => 1}) };
like($exception, qr/Can't open $metadata_path:/, "error when metadata missing and explicitly required");

$metadata = parse($opt, {required => 0});
is_deeply($metadata, {}, "metadata empty when missing and not required");

is(metaSampleName("filename_a", $opt), "sample", "default name for missing metadata (single-sample)");

is(HMF::Pipeline::Functions::Metadata::stripPath("/a/b/c.vcf", "b"), "c.vcf", "strip path to irrelevant part");
is(HMF::Pipeline::Functions::Metadata::stripPath("/a/b/c.vcf", "a"), "b/c.vcf", "strip path includes sub-directories");
is(HMF::Pipeline::Functions::Metadata::stripPath("/a/b/c.vcf", "t"), "/a/b/c.vcf", "strip path failure returns original path");

done_testing();
