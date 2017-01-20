#!/usr/bin/env perl

use discipline;

use version;
use Test::More;

use HMF::Pipeline;

# long story short: should always be developing the next version
# $VERSION == (current tag before "-") is allowed if exactly on current tag
# this is so that releases and release candidates pass tests
# $VERSION is ideally updated immediately after release to keep testes passing (but can be done whenever)

# once into the release cycle this test will begin to fail for *untagged* commits
# i.e. if you are releasing v1.11:
# it should pass on v1.11-rc.1, v1.11-rc.2, v1.11 etc.
# it should fail on commits between release candidates until you tag them
# it should fail on the next commit after v1.11 if this does not bump $VERSION

my $git_version = qx/git describe --tags/;
my $git_version_abbrev = qx/git describe --tags --abbrev=0/;
my ($git_major_version, $git_rest_version) = split "-", $git_version, 2;
my $git_parsed_version = version->parse($git_major_version);
my $package_version = version->parse($HMF::Pipeline::VERSION);
ok($package_version == $git_parsed_version && $git_version eq $git_version_abbrev || $package_version > $git_parsed_version, "package version up to date")
    or diag "$HMF::Pipeline::VERSION should be equal to ahead of $git_version";


done_testing();
