package HMF::Pipeline::Cobalt;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata qw(linkExtraArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING COBALT ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR}, cobalt => "cobalt");

    my @baf_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);

        say "COBALT Sample:" . $sample . " BAM:" . $sample_bam;
    }


    return;
}

1;