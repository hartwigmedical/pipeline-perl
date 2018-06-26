package HMF::Pipeline::PostStats;

use FindBin::libs;
use discipline;

use File::Spec::Functions qw(:ALL);

use HMF::Pipeline::Functions::Config qw(createDirs sampleControlBamsAndJobs sampleBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING POSTSTATS ###";

    my %snp_check_dirs = map { $_ => catfile($_, "snpcheck") } keys %{$opt->{SAMPLES}};
    my $dirs = createDirs($opt->{OUTPUT_DIR}, out => "QCStats", %snp_check_dirs);
    # SABR: weird relative path requirement

    my @designs;
    @designs = split '\t', $opt->{SNPCHECK_DESIGNS} if $opt->{SNPCHECK_DESIGNS};

    my $done_file = checkReportedDoneFile("PostStats", undef, $dirs, $opt) or return;

    my $qsub = qsubTemplate($opt, "POSTSTATS");
    my ($all_samples, $running_jobs) = sampleBamsAndJobs($opt);

    my $stats_job_id = fromTemplate(
        "PostStats",
        undef,
        0,
        $qsub,
        $running_jobs,
        $dirs,
        $opt,
        sample_bams => $all_samples,
        designs => \@designs,
    );

    my $check_job_id = fromTemplate(
        "PostStatsCheck",
        undef,
        0,
        $qsub,
        [$stats_job_id],
        $dirs,
        $opt,
        sample_bams => $all_samples,
        designs => \@designs,
    );
    my $done_job_id = markDone($done_file, [ $stats_job_id, $check_job_id ], $dirs, $opt);

    $opt->{RUNNING_JOBS}->{poststats} = [$done_job_id];

    return;
}

1;