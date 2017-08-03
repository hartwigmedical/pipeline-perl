package HMF::Pipeline::SomaticVariants;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs addSubDir getChromosomes sampleControlBamsAndJobs sampleBamAndJobs);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Job::Vcf qw(concat);
use HMF::Pipeline::Metadata;
use HMF::Pipeline::Sge qw(qsubTemplate qsubJava);
use HMF::Pipeline::BaseRecalibration qw(runRecalibrationOnSample);
use List::MoreUtils qw(uniq);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SOMATIC VARIANT CALLERS ###";

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "somaticVariants", $joint_name));

    my ($recalibrated_ref_bam, $recal_ref_jobs) = checkRecalibratedSample($ref_sample, $ref_bam_path, $opt);
    my ($recalibrated_tumor_bam, $recal_tumor_jobs) = checkRecalibratedSample($tumor_sample, $tumor_bam_path, $opt);
    $running_jobs = [ uniq @{$running_jobs}, @{$recal_ref_jobs}, @{$recal_tumor_jobs} ];

    say "\nRunning somatic callers on:";
    say "$joint_name \t $recalibrated_ref_bam \t $recalibrated_tumor_bam";

    my $done_file = checkReportedDoneFile($joint_name, undef, $dirs, $opt) or return;

    my ($job_id, $vcf) = runStrelka($ref_sample, $tumor_sample, $recalibrated_ref_bam, $recalibrated_tumor_bam, $joint_name, $running_jobs, $dirs, $opt);

    my $merge_job_ids = mergeSomatics($tumor_sample, $joint_name, $job_id, $vcf, $dirs, $opt);
    $job_id = markDone($done_file, [ $job_id, @{$merge_job_ids} ], $dirs, $opt);
    $opt->{RUNNING_JOBS}->{somvar} = [$job_id];

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

sub mergeSomatics {
    my ($tumor_sample, $joint_name, $strelka_job_id, $strelka_vcf, $dirs, $opt) = @_;

    say "\n### SCHEDULING MERGE SOMATIC VCFS ###";

    my @job_ids;
    my $qsub = qsubJava($opt, "SOMVARMERGE");
    my $output_vcf = $strelka_vcf;
    my $job_id;

    my $pre_annotate_vcf = $output_vcf;
    if ($opt->{SOMVAR_ANNOTATE} eq "yes") {
        (my $basename = $output_vcf) =~ s/\.vcf$//;
        $output_vcf = "${basename}_annotated.vcf";

        $job_id = fromTemplate(
            "SomaticAnnotation",
            undef,
            0,
            $qsub,
            [$strelka_job_id],
            $dirs,
            $opt,
            basename => $basename,
            final_vcf => $output_vcf,
        );
        push @job_ids, $job_id;
    }

    my $final_vcf = catfile($dirs->{out}, "${joint_name}_melted.vcf");
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

sub runStrelka {
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt) = @_;

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

1;
