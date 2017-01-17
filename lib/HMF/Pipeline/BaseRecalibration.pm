package HMF::Pipeline::BaseRecalibration;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Job::Bam;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING BASERECALIBRATION ###";
    say "Running base recalibration for the following BAM-files:";

    my $known_files = "";
    $known_files = join " ", map { "-knownSites $_" } split '\t', $opt->{BASERECALIBRATION_KNOWN} if $opt->{BASERECALIBRATION_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        HMF::Pipeline::Job::Bam::operationWithSliceChecks("BaseRecalibration", $sample, $known_files, "recalibrated", "recal", $opt);
    }
    return;
}

1;
