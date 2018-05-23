package HMF::Pipeline::PostStats;

use FindBin::libs;
use discipline;

use File::Spec::Functions qw(:ALL);

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Metadata qw(linkArtefact linkExtraArtefact);

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

    my ($sample_bams, $running_jobs) = sampleBamsAndJobs($opt);
    my $qsub = qsubTemplate($opt, "POSTSTATS");
    my $done_file = checkReportedDoneFile("PostStats", undef, $dirs, $opt) or return;
    my $stats_job_id = fromTemplate(
        "PostStats",
        undef,
        0,
        $qsub,
        $running_jobs,
        $dirs,
        $opt,
        sample_bams => $sample_bams,
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
        sample_bams => $sample_bams,
        designs => \@designs,
    );
    my $job_id = markDone($done_file, [ $stats_job_id, $check_job_id ], $dirs, $opt);

    $opt->{RUNNING_JOBS}->{poststats} = [$job_id];
    linkArtefacts(\@designs, $dirs, $opt);

    return;
}

sub linkArtefacts {
    my ($designs, $dirs, $opt) = @_;

    # SABR: dependent on implicit bamMetrics naming
    my $metrics_path = catfile($dirs->{out}, "$opt->{RUN_NAME}.bamMetrics.pdf");
    linkArtefact($metrics_path, "qc", $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        # SABR: dependent on implicit SNPcheck naming (easier to fix)
        foreach my $design (@{$designs}) {
            linkExtraArtefact(catfile($dirs->{$sample}, "${sample}_${design}.vcf"), $opt);
        }
    }
    return;
}

1;