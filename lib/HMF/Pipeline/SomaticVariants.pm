package HMF::Pipeline::SomaticVariants;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs addSubDir sampleControlBamsAndJobs);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Metadata;
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::BaseRecalibration qw(runRecalibrationOnSample);
use List::MoreUtils qw(uniq);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SOMATIC VARIANT CALLING ###";
    $opt->{RUNNING_JOBS}->{somvar} = [];

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "somaticVariants", $joint_name));
    my $final_vcf = catfile($dirs->{out}, "${joint_name}_post_processed.vcf");
    $opt->{SOMVAR_VCF_FILE} = $final_vcf;

    my ($recalibrated_ref_bam, $recal_ref_jobs) = checkRecalibratedSample($ref_sample, $ref_bam_path, $opt);
    my ($recalibrated_tumor_bam, $recal_tumor_jobs) = checkRecalibratedSample($tumor_sample, $tumor_bam_path, $opt);
    $running_jobs = [ uniq @{$running_jobs}, @{$recal_ref_jobs}, @{$recal_tumor_jobs} ];

    say "\nRunning somatic calling on:";
    say "$joint_name \t $recalibrated_ref_bam \t $recalibrated_tumor_bam";
    my $done_file = checkReportedDoneFile("Somatic_$joint_name", undef, $dirs, $opt) or return;

    my ($job_id, $vcf) = runStrelka($tumor_sample, $recalibrated_ref_bam, $recalibrated_tumor_bam, $joint_name, $running_jobs, $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{somvar}}, $job_id;

    my $post_process_job_ids = postProcessStrelka($final_vcf, $job_id, $vcf, $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{somvar}}, @{$post_process_job_ids};

    $job_id = markDone($done_file, [ $job_id, @{$post_process_job_ids} ], $dirs, $opt);
    push @{$opt->{RUNNING_JOBS}->{somvar}}, $job_id;
    return;
}

sub checkRecalibratedSample {
    my ($sample, $sample_bam_path, $opt) = @_;

    if (index($sample_bam_path, "recalibrated.bam") == -1) {
        say "\nMissing recalibrated file for sample: $sample";
        my ($recalibrated_bam, $recalibration_jobs) = HMF::Pipeline::BaseRecalibration::runRecalibrationOnSample($sample, $opt);
        return ($recalibrated_bam, $recalibration_jobs);
    }
    return ($sample_bam_path, []);
}

sub runStrelka {
    my ($tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING STRELKA ###";

    $dirs->{strelka}->{out} = addSubDir($dirs, "strelka");
    my $final_vcf = catfile($dirs->{strelka}->{out}, "passed.somatic.merged.vcf");

    my $job_id = fromTemplate(
        "Strelka",
        undef,
        1,
        qsubJava($opt, "STRELKA"),
        $running_jobs,
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        joint_name => $joint_name,
        ref_bam_path => $ref_bam_path,
        tumor_bam_path => $tumor_bam_path,
        final_vcf => $final_vcf,
    );

    return ($job_id, $final_vcf);
}

sub postProcessStrelka {
    my ($final_vcf, $strelka_job_id, $strelka_vcf, $dirs, $opt) = @_;

    say "\n### SCHEDULING STRELKA POST PROCESS ###";

    my @job_ids;
    my $qsub = qsubJava($opt, "SOMVARPOSTPROCESS");
    my $output_vcf = $strelka_vcf;
    my $job_id = $strelka_job_id;

    my $pre_annotate_vcf = $output_vcf;
    if ($opt->{SOMVAR_ANNOTATE} eq "yes") {
        (my $basename = $output_vcf) =~ s/\.vcf$//;
        $output_vcf = "${basename}_annotated.vcf";

        $job_id = fromTemplate(
            "SomaticAnnotation",
            undef,
            0,
            $qsub,
            [$job_id],
            $dirs,
            $opt,
            basename => $basename,
            final_vcf => $output_vcf,
        );
        push @job_ids, $job_id;
    }

    $job_id = fromTemplate(
        "SomaticPONAnnotation",
        undef,
        0,
        $qsub,
        [$job_id],
        $dirs,
        $opt,
        pre_annotate_vcf => $pre_annotate_vcf,
        input_vcf => $output_vcf,
        output_vcf => $final_vcf,
    );
    push @job_ids, $job_id;

    HMF::Pipeline::Metadata::linkVcfArtefacts($final_vcf, "somatic", $opt) if $job_id;

    return \@job_ids;
}

1;
