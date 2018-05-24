#!/usr/bin/env perl

use FindBin::libs;
use discipline;

use HMF::Pipeline::Functions::Config;
use HMF::Pipeline;

my $opt = {};
die usage() if @ARGV == 0;
HMF::Pipeline::Functions::Config::parse($ARGV[0], $opt);
HMF::Pipeline::Functions::Config::validate($opt);
HMF::Pipeline::Functions::Config::createDirs($opt->{OUTPUT_DIR});
HMF::Pipeline::Functions::Config::setupLogging($opt->{OUTPUT_DIR});
HMF::Pipeline::lockRun($opt->{OUTPUT_DIR});
HMF::Pipeline::Functions::Config::addSamples($opt);
HMF::Pipeline::Functions::Config::recordGitVersion($opt);
HMF::Pipeline::Functions::Config::copyConfigAndScripts($opt);
HMF::Pipeline::run($opt);

sub usage {
    warn <<"END";
Usage: $0 configurationFile.conf
END
    exit;
}
