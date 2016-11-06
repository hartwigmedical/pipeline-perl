#!/usr/bin/env perl

use strict;
use warnings;

use version;
use Test::More;

use HMF::Pipeline;


my $git_version = qx/git describe --tags/;
my ($git_major_version, $git_rest_version) = split "-", $git_version, 2;
my $git_parsed_version = version->parse($git_major_version);
my $package_version = version->parse($HMF::Pipeline::VERSION);
ok($package_version == $git_parsed_version && !$git_rest_version || $package_version > $git_parsed_version, "package version up to date")
    or diag "$HMF::Pipeline::VERSION should be equal to ahead of $git_version";


done_testing();
