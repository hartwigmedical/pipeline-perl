package HMF::Pipeline::StructuralVariants;

use FindBin::libs;
use discipline;

use File::Spec::Functions;
use File::Basename;
use Sort::Key::Natural qw(mkkey_natural);

use HMF::Pipeline::Functions::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkVcfArtefacts);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING STRUCTURAL VARIANT CALLING ###";

    $opt->{RUNNING_JOBS}->{'sv'} = [];
    if ($opt->{MANTA} eq "yes") {
        my $manta_jobs = runManta($opt);
        push @{$opt->{RUNNING_JOBS}->{'sv'}}, @{$manta_jobs};
    }

    if ($opt->{GRIDSS} eq "yes") {
        my $gridss_jobs = runGridss($opt);
        if ($gridss_jobs) {
            push @{$opt->{RUNNING_JOBS}->{'sv'}}, @{$gridss_jobs};
        }
    }

    return;
}

sub runGridss {
    my ($opt) = @_;

    say "\n### SCHEDULING GRIDSS ###";

    my @gridss_jobs;
    my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, undef) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "gridss", $joint_name));

    my $done_file = checkReportedDoneFile("Gridss_$joint_name", undef, $dirs, $opt) or return;

    # KODU: Poststats depends on the BAM creation, so fine to depend on the poststats job.
    my ($ref_pre_process_job_id, undef) =
        runGridssPreProcess($dirs, $ref_sample, $ref_sample_bam, $opt->{REF_INSERT_SIZE_METRICS}, $opt->{RUNNING_JOBS}->{poststats}, $opt);
    push @gridss_jobs, $ref_pre_process_job_id;
    my ($tumor_pre_process_job_id, undef) =
        runGridssPreProcess($dirs, $tumor_sample, $tumor_sample_bam, $opt->{TUMOR_INSERT_SIZE_METRICS}, $opt->{RUNNING_JOBS}->{poststats}, $opt);
    push @gridss_jobs, $tumor_pre_process_job_id;

    my $done_job_id = markDone($done_file, [ $ref_pre_process_job_id, $tumor_pre_process_job_id ], $dirs, $opt);
    push @gridss_jobs, $done_job_id;

    return \@gridss_jobs;
}

sub runGridssPreProcess {
    my ($dirs, $sample, $sample_bam, $insert_size_metrics, $dependent_jobs, $opt) = @_;
    my $working_dir = catfile($dirs->{out}, join "", basename($sample_bam), ".gridss.working");
    my $pre_process_bam = catfile($working_dir, join "", basename($sample_bam), ".sv.bam");

    my $job_id = fromTemplate(
        "GridssPreProcess",
        $sample,
        1,
        qsubTemplate($opt, "GRIDSS"),
        $dependent_jobs,
        $dirs,
        $opt,
        sample => basename($sample_bam),
        sample_bam => $sample_bam,
        insert_size_metrics => $insert_size_metrics,
        working_dir => $working_dir,
        pre_process_bam => $pre_process_bam,
    );

    return ($job_id, $pre_process_bam);
}

sub runManta {
    my ($opt) = @_;

    say "\n### SCHEDULING MANTA ###";

    my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);

    my @manta_jobs;
    my $job_id = runMantaJob($tumor_sample_bam, $ref_sample_bam, $joint_name, $running_jobs, $opt);
    push @manta_jobs, $job_id;

    $job_id = runBreakpointInspector($tumor_sample, $tumor_sample_bam, $ref_sample, $ref_sample_bam, $joint_name, $job_id, $opt);
    push @manta_jobs, $job_id;
    return \@manta_jobs;
}

sub runMantaJob {
    my ($sample_bam, $control_bam, $joint_name, $running_jobs, $opt) = @_;

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "manta", $joint_name));

    my $job_id = fromTemplate(
        "Manta",
        undef,
        1,
        qsubTemplate($opt, "MANTA"),
        $running_jobs,
        $dirs,
        $opt,
        sample_bam => $sample_bam,
        control_bam => $control_bam,
        joint_name => $joint_name,
    );

    return $job_id;
}

sub runBreakpointInspector {
    my ($sample, $sample_bam, $control, $control_bam, $joint_name, $manta_job_id, $opt) = @_;

    my $manta_vcf = catfile($opt->{OUTPUT_DIR}, "structuralVariants", "manta", $joint_name, "results", "variants", "somaticSV.vcf.gz");

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "bpi", $joint_name));
    $opt->{BPI_VCF_FILE} = catfile($dirs->{out}, "${joint_name}_somaticSV_bpi.vcf");

    my $job_id = fromTemplate(
        "BreakpointInspector",
        undef,
        1,
        qsubTemplate($opt, "BPI"),
        [$manta_job_id],
        $dirs,
        $opt,
        sample => $sample,
        control => $control,
        sample_bam => $sample_bam,
        control_bam => $control_bam,
        joint_name => $joint_name,
        input_vcf => $manta_vcf,
    );

    $opt->{BPI_VCF_FILE} = join "", $opt->{BPI_VCF_FILE}, ".gz";
    linkVcfArtefacts($opt->{BPI_VCF_FILE}, 'structural_variant', $opt);

    return $job_id;
}

1;
