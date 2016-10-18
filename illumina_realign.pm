package illumina_realign;

use 5.16.0;
use strict;
use warnings;

use File::Basename;
use File::Spec::Functions;

use FindBin;
use lib "$FindBin::Bin";

use illumina_jobs qw(bamOperationWithSliceChecks);


sub runRealignment {
    my ($opt) = @_;

    say "Running single sample indel realignment for the following BAM-files:";

    my $known_files;
    $known_files = join " ", map { "-known $_" } split '\t', $opt->{REALIGNMENT_KNOWN} if $opt->{REALIGNMENT_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        bamOperationWithSliceChecks("Realignment", $sample, $known_files, "realigned", "realign", $opt);
    }
    return;
}

1;
