#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Test::More;

foreach my $template (glob("templates/*.sh.tt")) {
    open my $fh, "<", $template;
    my @lines = grep { /INCLUDE Logging.tt/ } <$fh>;
    close $fh;
    is(@lines, 1, "precisely one include of logging functions");
    my $job_name = fileparse($template, qr/\.sh\.tt$/);
    like($lines[0], qr/\sjob_name="$job_name"\s/, "logged job name matches actual job name");
}

done_testing();
