package HMF::Pipeline::Cobalt;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING COBALT ANALYSIS ###";

    my @cobalt_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);
        my $dirs = createDirs($opt->{OUTPUT_DIR}, cobalt => "cobalt");

        my $job_id = fromTemplate(
            "Cobalt",
            $sample,
            1,
            qsubJava($opt, "COBALT"),
            $running_jobs,
            $dirs,
            $opt,
            sample => $sample,
            sample_bam => $sample_bam,
        );
        next unless $job_id;

        push @cobalt_jobs, $job_id;
    }
    $opt->{RUNNING_JOBS}->{'cobalt'} = \@cobalt_jobs;
    return;
}

1;