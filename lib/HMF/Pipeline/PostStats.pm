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

    # KODU: We need the insert size metrics when running gridss. Their naming comes out of poststats and is dependent on the mode we run in.
    my ($ref_sample, $tumor_sample, undef, undef, undef, undef) = sampleControlBamsAndJobs($opt);
    my $suffix = "_MultipleMetrics.txt.insert_size_metrics";
    my $ref_sample_name;
    my $tumor_sample_name;
    if ($opt->{BAM}) {
        $ref_sample_name = $ref_sample;
        $tumor_sample_name = $tumor_sample;
    } else { # KODU: Run from FASTQ
        $ref_sample_name = join "", $ref_sample, "_dedup";
        $tumor_sample_name = join "", $tumor_sample, "_dedup";
    }

    $opt->{REF_INSERT_SIZE_METRICS} = catfile($opt->{OUTPUT_DIR}, "QCStats", $ref_sample_name, join "", $ref_sample_name, $suffix);
    $opt->{TUMOR_INSERT_SIZE_METRICS} = catfile($opt->{OUTPUT_DIR}, "QCStats", $tumor_sample_name, join "", $tumor_sample_name, $suffix);

    say join "", "\nRef insert size: ", $opt->{REF_INSERT_SIZE_METRICS};
    say join "", "\nTumor insert size: ", $opt->{TUMOR_INSERT_SIZE_METRICS};

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
    my $job_id = markDone($done_file, [ $stats_job_id, $check_job_id ], $dirs, $opt);

    $opt->{RUNNING_JOBS}->{poststats} = [$job_id];

    return;
}

1;