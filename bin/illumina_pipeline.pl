#!/usr/bin/env perl

use FindBin::libs;
use discipline;

use illumina_config qw(readConfig checkConfig setupLogging recordGitVersion copyConfigAndScripts);
use HMF::Pipeline qw(getSamples createOutputDirs lockRun runPipeline);


my $opt = {};
die usage() if @ARGV == 0;
readConfig($ARGV[0], $opt);
checkConfig($opt);
getSamples($opt);
createOutputDirs($opt->{OUTPUT_DIR}, $opt->{SAMPLES});
setupLogging($opt->{OUTPUT_DIR});
lockRun($opt->{OUTPUT_DIR});
recordGitVersion($opt);
copyConfigAndScripts($opt);
runPipeline($opt);


sub usage {
    warn <<"END";
    Usage: perl illumina_pipeline.pl configurationFile.conf
END
    exit;
}
