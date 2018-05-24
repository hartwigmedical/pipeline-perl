package HMF::Pipeline::Realignment;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Bam;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING INDEL REALIGNMENT ###";
    say "Running single sample indel realignment for the following BAM-files:";

    my $known_files = "";
    $known_files = join " ", map { "-known $_" } split '\t', $opt->{REALIGNMENT_KNOWN} if $opt->{REALIGNMENT_KNOWN};

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        HMF::Pipeline::Job::Bam::operationWithSliceChecks("Realignment", $sample, $known_files, "realigned", "realign", $opt);
    }
    return;
}

1;