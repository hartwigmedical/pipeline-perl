package HMF::Pipeline::StructuralVariants;

use FindBin::libs;
use discipline;

use File::Spec::Functions;
use Sort::Key::Natural qw(mkkey_natural);

use HMF::Pipeline::Functions::Config qw(createDirs sampleBamAndJobs sampleBamsAndJobs sampleControlBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate checkReportedDoneFile);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkArtefact linkExtraArtefact);

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
    return;
}

sub runManta {
    my ($opt) = @_;

    say "\n### SCHEDULING MANTA ###";

    my @manta_jobs;
    my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    say "\n$joint_name \t $ref_sample_bam \t $tumor_sample_bam";
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

    linkExtraArtefact(catfile($dirs->{out}, "results", "variants", "diploidSV.vcf.gz"), $opt);
    linkExtraArtefact(catfile($dirs->{out}, "results", "variants", "diploidSV.vcf.gz.tbi"), $opt);
    if (defined $control_bam) {
        linkExtraArtefact(catfile($dirs->{out}, "results", "variants", "somaticSV.vcf.gz"), $opt);
        linkExtraArtefact(catfile($dirs->{out}, "results", "variants", "somaticSV.vcf.gz.tbi"), $opt);
    }

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

    linkArtefact($opt->{BPI_VCF_FILE}, 'somatic_sv', $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${control}_sliced.bam"), $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${sample}_sliced.bam"), $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${joint_name}_bpi_stats.tsv"), $opt);

    return $job_id;
}

1;
