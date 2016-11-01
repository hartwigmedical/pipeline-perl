package illumina_baseRecal;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;

use illumina_jobs qw(bamOperationWithSliceChecks);

use parent qw(Exporter);
our @EXPORT_OK = qw(runBaseRecalibration);


sub runBaseRecalibration {
    my ($opt) = @_;

    say "\n### SCHEDULING BASERECALIBRATION ###";
    say "Running base recalibration for the following BAM-files:";

    my $known_files = "";
    $known_files = join " ", map { "-knownSites $_" } split '\t', $opt->{BASERECALIBRATION_KNOWN} if $opt->{BASERECALIBRATION_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        bamOperationWithSliceChecks("BaseRecalibration", $sample, $known_files, "recalibrated", "recal", $opt);
    }
    return;
}

1;
