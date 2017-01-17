#!/usr/bin/env perl

use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec::Functions;
use File::Temp;
use Test::Fatal;
use Test::More;

use HMF::Pipeline::Metadata qw(parse metaSampleName sampleControlNames linkArtefact linkExtraArtefact linkVcfArtefacts linkBamArtefacts writeLinks);


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

my $temp_dir = File::Temp->newdir();
my $opt = {OUTPUT_DIR => $temp_dir};

my $metadata_path = catfile($temp_dir, "metadata");
HMF::Pipeline::Metadata::writeJson(
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
        type_a_bam => catfile($temp_dir, "sample_a", "mapping", "filename_a.bam"),
        type_a_bai => catfile($temp_dir, "sample_a", "mapping", "filename_a.bam.bai"),
        type_b_bam => catfile($temp_dir, "sample_b", "mapping", "filename_b.bam"),
        type_b_bai => catfile($temp_dir, "sample_b", "mapping", "filename_b.bam.bai"),
    },
    "bam filename links stored"
);

delete $opt->{LINKS};
linkArtefact("filename_a", "artefact_a", $opt);
linkArtefact("filename_b", "artefact_b", $opt);
ok(exists $opt->{LINKS}, "links stored in \$opt");
is_deeply($opt->{LINKS}, {artefact_a => "filename_a", artefact_b => "filename_b"}, "artefact filename links stored");

linkExtraArtefact(catfile($temp_dir, "filename_a"), $opt);
linkExtraArtefact(catfile($temp_dir, "filename_b"), $opt);
ok(exists $opt->{EXTRAS}, "extra links stored in \$opt");
is_deeply($opt->{EXTRAS}, [ "filename_a", "filename_b" ], "extra artefact filenames stored");

make_path(catfile($temp_dir, "logs"));
writeLinks($opt);
my $links_path = catfile($temp_dir, "logs", "links.json");
my $links = HMF::Pipeline::Metadata::readJson($links_path);
is_deeply(
    $links, {
        artefact_a => "filename_a",
        artefact_b => "filename_b",
    },
    "links written to file"
);

delete $opt->{LINKS};
linkVcfArtefacts("filename", "artefact", $opt);
ok(exists $opt->{LINKS}, "extra links stored in \$opt");
is_deeply($opt->{LINKS}, {artefact_vcf => "filename", artefact_vcf_index => "filename.idx"}, "vcf artefact links stored");


HMF::Pipeline::Metadata::writeJson(
    $metadata_path, {
        ref_sample => "ref",
        tumor_sample => "tumor",
    }
);
my $names = [ sampleControlNames($opt) ];
is_deeply($names, [ "ref", "tumor", "ref_tumor" ], "gets sample/control/joint names");

my $exception;
HMF::Pipeline::Metadata::writeJson(
    $metadata_path, {
        tumor_sample => "tumor",
    }
);
$exception = exception { my $names = [ sampleControlNames($opt) ] };
like($exception, qr/metadata missing ref_sample/, "detects missing ref sample");

HMF::Pipeline::Metadata::writeJson(
    $metadata_path, {
        ref_sample => "ref",
    }
);
$exception = exception { my $names = [ sampleControlNames($opt) ] };
like($exception, qr/metadata missing tumor_sample/, "detects missing tumor sample");

$temp_dir->DESTROY();
$exception = exception { $metadata = HMF::Pipeline::Metadata::readJson($links_path) };
like($exception, qr/Can't open $links_path:/, "fails to open missing file");
$exception = exception { HMF::Pipeline::Metadata::writeJson($links_path, {}) };
like($exception, qr/Can't open $links_path:/, "fails to write to non-existent directory");

$exception = exception { $metadata = parse($opt) };
like($exception, qr/Can't open $metadata_path:/, "error when metadata missing");

$exception = exception { $metadata = parse($opt, {required => 1}) };
like($exception, qr/Can't open $metadata_path:/, "error when metadata missing and explicitly required");

$metadata = parse($opt, {required => 0});
is_deeply($metadata, {}, "metadata empty when missing and not required");

is(metaSampleName("filename_a", $opt), "sample", "default name for missing metadata (single-sample)");

done_testing();
