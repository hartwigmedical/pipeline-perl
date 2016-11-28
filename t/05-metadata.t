#!/usr/bin/env perl

use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec::Functions;
use File::Temp;
use Test::Fatal;
use Test::More;

use HMF::Pipeline::Metadata qw(parse metaSampleName linkArtefact linkExtraArtefact writeLinks);


## no critic (Subroutines::ProhibitCallsToUnexportedSubs)

my $temp_dir = File::Temp->newdir();
my $opt = {OUTPUT_DIR => $temp_dir};

my $metadata_path = catfile($temp_dir, "metadata");
HMF::Pipeline::Metadata::writeJson(
    $metadata_path, {
        test_sample_a => "filename_a",
        test_sample_b => "filename_b",
    }
);

my $metadata = parse($opt);
is_deeply(
    $metadata, {
        test_sample_a => "filename_a",
        test_sample_b => "filename_b",
    },
    "reads metadata"
);

is(metaSampleName("filename_a", $opt), "test_sample_a", "names sample file from metadata");
is(metaSampleName("filename_b", $opt), "test_sample_b", "names another sample file from metadata");
is(metaSampleName("filename_c", $opt), "sample", "default name for missing sample file");

linkArtefact("filename_a", "artefact_a", $opt);
linkArtefact("filename_b", "artefact_b", $opt);
ok(exists $opt->{LINKS}, "links stored in \$opt");
is($opt->{LINKS}->{artefact_a}, "filename_a", "artefact filename stored");
is($opt->{LINKS}->{artefact_b}, "filename_b", "another artefact filename stored");

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

my $exception;
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
