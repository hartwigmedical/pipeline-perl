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
    my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, $running_sample_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "gridss", $joint_name));

    my $done_file = checkReportedDoneFile("Gridss_$joint_name", undef, $dirs, $opt) or return;

    # KODU: GRIDSS requires the insert size metrics output from poststats, so should wait on poststats to finish.
    my $dependent_jobs = [ @{$running_sample_jobs}, @{$opt->{RUNNING_JOBS}->{poststats}} ];
    say $dependent_jobs;
    my ($ref_pre_process_job_id) =
        runGridssPreProcess($dirs, $ref_sample, $ref_sample_bam, $opt->{REF_INSERT_SIZE_METRICS}, $dependent_jobs, $opt);
    push @gridss_jobs, $ref_pre_process_job_id;
    my ($tumor_pre_process_job_id) =
        runGridssPreProcess($dirs, $tumor_sample, $tumor_sample_bam, $opt->{TUMOR_INSERT_SIZE_METRICS}, $dependent_jobs, $opt);
    push @gridss_jobs, $tumor_pre_process_job_id;

    my ($assemble_job_id, $assembly_bam) = runGridssAssemble($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, \@gridss_jobs, $opt);
    push @gridss_jobs, $assemble_job_id;

    my ($calling_job_id, $gridss_raw_vcf) = runGridssCalling($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, $assembly_bam, \@gridss_jobs, $opt);
    push @gridss_jobs, $calling_job_id;

    my ($annotation_job_id) = runGridssAnnotation($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, $assembly_bam, $gridss_raw_vcf, \@gridss_jobs, $opt);
    push @gridss_jobs, $annotation_job_id;

    my $done_job_id = markDone($done_file, \@gridss_jobs, $dirs, $opt);
    push @gridss_jobs, $done_job_id;

    return \@gridss_jobs;
}

sub runGridssPreProcess {
    my ($dirs, $sample, $sample_bam, $insert_size_metrics, $dependent_jobs, $opt) = @_;
    my $working_dir = catfile($dirs->{out}, join "", basename($sample_bam), ".gridss.working");
    my $sample_sv_bam = catfile($working_dir, join "", basename($sample_bam), ".sv.bam");

    my $job_id = fromTemplate(
        "GridssPreProcess",
        $sample,
        1,
        qsubTemplate($opt, "GRIDSS_PREPROCESS"),
        $dependent_jobs,
        $dirs,
        $opt,
        sample => basename($sample_bam),
        sample_bam => $sample_bam,
        insert_size_metrics => $insert_size_metrics,
        working_dir => $working_dir,
        sv_bam => $sample_sv_bam,
    );

    return ($job_id);
}

sub runGridssAssemble {
    my ($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, $dependent_jobs, $opt) = @_;

    my $assembly_bam_name = join "", $joint_name, ".assembly.bam";
    my $working_dir_name = join "", $assembly_bam_name, ".gridss.working";

    my $assembly_bam = catfile($dirs->{out}, $assembly_bam_name);
    my $metrics_output_dir = catfile($dirs->{out}, $working_dir_name);
    my $metrics_output = catfile($metrics_output_dir, $assembly_bam_name);
    my $assembly_sv_bam = catfile($dirs->{out}, $working_dir_name, join "", $assembly_bam_name, ".sv.bam");

    my $job_id = fromTemplate(
        "GridssAssemble",
        undef,
        1,
        qsubTemplate($opt, "GRIDSS_ASSEMBLE"),
        $dependent_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        joint_name => $joint_name,
        assembly_bam => $assembly_bam,
        metrics_output_dir => $metrics_output_dir,
        metrics_output => $metrics_output,
        sv_bam => $assembly_sv_bam
    );

    return ($job_id, $assembly_bam);
}

sub runGridssCalling {
    my ($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, $assembly_bam, $dependent_jobs, $opt) = @_;

    my $gridss_raw_vcf = catfile($dirs->{out}, join "", $joint_name, ".raw.gridss.vcf");

    my $job_id = fromTemplate(
        "GridssCalling",
        undef,
        1,
        qsubTemplate($opt, "GRIDSS_CALLING"),
        $dependent_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        joint_name => $joint_name,
        assembly_bam => $assembly_bam,
        gridss_raw_vcf => $gridss_raw_vcf
    );

    return ($job_id, $gridss_raw_vcf);
}

sub runGridssAnnotation {
    my ($dirs, $ref_sample_bam, $tumor_sample_bam, $joint_name, $assembly_bam, $gridss_raw_vcf, $dependent_jobs, $opt) = @_;

    my $gridss_annotated_vcf = catfile($dirs->{out}, join "", $joint_name, ".gridss.vcf");

    my $job_id = fromTemplate(
        "GridssAnnotation",
        undef,
        1,
        qsubTemplate($opt, "GRIDSS_ANNOTATE"),
        $dependent_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        joint_name => $joint_name,
        assembly_bam => $assembly_bam,
        gridss_raw_vcf => $gridss_raw_vcf,
        gridss_annotated_vcf => $gridss_annotated_vcf
    );

    return ($job_id);
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
    my ($tumor_sample_bam, $ref_sample_bam, $joint_name, $running_jobs, $opt) = @_;

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "manta", $joint_name));

    my $job_id = fromTemplate(
        "Manta",
        undef,
        1,
        qsubTemplate($opt, "MANTA"),
        $running_jobs,
        $dirs,
        $opt,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        joint_name => $joint_name,
    );

    return $job_id;
}

sub runBreakpointInspector {
    my ($sample, $sample_bam, $control, $control_bam, $joint_name, $manta_job_id, $opt) = @_;

    my $manta_vcf = catfile($opt->{OUTPUT_DIR}, "structuralVariants", "manta", $joint_name, "results", "variants", "somaticSV.vcf.gz");

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "bpi", $joint_name));
    $opt->{STRUCTURAL_VARIANT_VCF} = catfile($dirs->{out}, "${joint_name}_somaticSV_bpi.vcf");

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

    $opt->{STRUCTURAL_VARIANT_VCF} = join "", $opt->{STRUCTURAL_VARIANT_VCF}, ".gz";
    linkVcfArtefacts($opt->{STRUCTURAL_VARIANT_VCF}, 'structural_variant', $opt);

    return $job_id;
}

1;
