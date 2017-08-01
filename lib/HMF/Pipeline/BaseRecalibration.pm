package HMF::Pipeline::BaseRecalibration;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Job::Bam;

use parent qw(Exporter);
our @EXPORT_OK = qw(run
    runRecalibrationOnSample
);


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

sub runRecalibrationOnSample {
    my ($sample, $opt) = @_;
    say "\n### SCHEDULING BASERECALIBRATION ###";
    say "Running base recalibration for the following BAM:";

    my $known_files = "";
    $known_files = join " ", map { "-knownSites $_" } split '\t', $opt->{BASERECALIBRATION_KNOWN} if $opt->{BASERECALIBRATION_KNOWN};

    my $sample_bam = $opt->{BAM_FILES}->{$sample};

    my ($recalibrated_bam, $job_id) = HMF::Pipeline::Job::Bam::bamOperationWithSliceChecks("BaseRecalibration", $sample, $sample_bam, $known_files, "recalibrated", "recal", $opt);
    return ($recalibrated_bam, $job_id);
}

1;
