package HMF::Pipeline::PostStats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions qw(:ALL);

use HMF::Pipeline::Config qw(createDirs sampleBamsAndJobs);
use HMF::Pipeline::Job qw(getId fromTemplate);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Metadata qw(linkArtefact linkExtraArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING POSTSTATS ###";

    my %snp_check_dirs = map { $_ => catfile($_, "snpcheck") } keys %{$opt->{SAMPLES}};
    my $dirs = createDirs($opt->{OUTPUT_DIR}, out => "QCStats", %snp_check_dirs);
    # weird relative path requirement
    my $exoncov_dir = catfile($dirs->{out}, "exoncov");
    $dirs->{exoncov} = abs2rel($exoncov_dir, $dirs->{tmp});

    my @designs;
    @designs = split '\t', $opt->{SNPCHECK_DESIGNS} if $opt->{SNPCHECK_DESIGNS};

    my ($sample_bams, $running_jobs) = sampleBamsAndJobs($opt);
    my $qsub = qsubTemplate($opt, "POSTSTATS");
    my $job_id = fromTemplate(
        "PostStats",
        undef,
        1,
        $qsub,
        $running_jobs,
        $dirs,
        $opt,
        sample_bams => $sample_bams,
        designs => \@designs,
    );

    linkArtefacts(\@designs, $exoncov_dir, $dirs, $opt);
    return if not $job_id;

    my $job_id_check = fromTemplate(
        "PostStatsCheck",
        undef,
        0,
        $qsub,
        [$job_id],
        $dirs,
        $opt,
        sample_bams => $sample_bams,
        designs => \@designs,
    );

    $opt->{RUNNING_JOBS}->{poststats} = [$job_id_check] if $job_id_check;
    return;
}

sub linkArtefacts {
    my ($designs, $exoncov_dir, $dirs, $opt) = @_;

    # dependent on implicit bamMetrics naming
    my $metrics_path = catfile($dirs->{out}, "$opt->{RUN_NAME}.bamMetrics.pdf");
    linkArtefact($metrics_path, "qc", $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        # dependent on implicit SNPcheck naming (easier to fix)
        foreach my $design (@{$designs}) {
            linkExtraArtefact(catfile($dirs->{$sample}, "${sample}_${design}.vcf"), $opt);
        }

        if ($opt->{EXONCALLCOV} eq "yes") {
            # dependent on implicit ExonCov naming
            linkExtraArtefact(catfile($exoncov_dir, "${sample}.html"), $opt);
            (my $coverage_name = $opt->{BAM_FILES}->{$sample}) =~ s/\.bam$/_exon_coverage.tsv/;
            linkExtraArtefact(catfile($exoncov_dir, "Exons", $coverage_name), $opt);
        }
    }
    return;
}

1;
