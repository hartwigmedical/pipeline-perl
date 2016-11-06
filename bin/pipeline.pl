#!/usr/bin/env perl

use FindBin::libs;
use discipline;

use HMF::Pipeline::Config;
use HMF::Pipeline;


my $opt = {};
die usage() if @ARGV == 0;
HMF::Pipeline::Config::parse($ARGV[0], $opt);
HMF::Pipeline::Config::validate($opt);
HMF::Pipeline::Config::createDirs($opt->{OUTPUT_DIR});
HMF::Pipeline::Config::setupLogging($opt->{OUTPUT_DIR});
HMF::Pipeline::lockRun($opt->{OUTPUT_DIR});
HMF::Pipeline::Config::addSamples($opt);
HMF::Pipeline::Config::recordGitVersion($opt);
HMF::Pipeline::Config::copyConfigAndScripts($opt);
HMF::Pipeline::run($opt);


sub usage {
    warn <<"END";
Usage: $0 configurationFile.conf
END
    exit;
}
