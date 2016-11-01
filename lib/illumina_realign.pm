package illumina_realign;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_jobs qw(bamOperationWithSliceChecks);

use parent qw(Exporter);
our @EXPORT_OK = qw(runRealignment);


sub runRealignment {
    my ($opt) = @_;

    say "\n### SCHEDULING INDELREALIGNMENT ###";
    say "Running single sample indel realignment for the following BAM-files:";

    my $known_files = "";
    $known_files = join " ", map { "-known $_" } split '\t', $opt->{REALIGNMENT_KNOWN} if $opt->{REALIGNMENT_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        bamOperationWithSliceChecks("Realignment", $sample, $known_files, "realigned", "realign", $opt);
    }
    return;
}

1;
