package HMF::Pipeline::SomaticVariants;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs addSubDir getChromosomes sampleControlBamsAndJobs);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Job::Vcf qw(concat);
use HMF::Pipeline::Metadata;
use HMF::Pipeline::Sge qw(qsubTemplate qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SOMATIC VARIANT CALLERS ###";

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "somaticVariants", $joint_name));

    if (index($tumor_bam_path, "recalibrated.bam") == -1) {
        HMF::Pipeline::BaseRecalibration::runRecalibrationOnSample($tumor_sample, $opt);
        my ($tumor_sample_bam, $tumor_sample_jobs) = sampleBamAndJobs($tumor_sample, $opt);
        $tumor_bam_path = $tumor_sample_bam;
        $running_jobs = [ uniq @{$running_jobs}, @{$tumor_sample_jobs} ];
    }

    say "\n$joint_name \t $ref_bam_path \t $tumor_bam_path";

    my $done_file = checkReportedDoneFile($joint_name, undef, $dirs, $opt) or return;

    my ($job_id, $vcf) = runStrelka($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt);

    my $merge_job_ids = mergeSomatics($tumor_sample, $joint_name, $job_id, $vcf, $dirs, $opt);
    $job_id = markDone($done_file, [ $job_id, @{$merge_job_ids} ], $dirs, $opt);
    $opt->{RUNNING_JOBS}->{somvar} = [$job_id];

    return;
}

sub mergeSomatics {
    my ($tumor_sample, $joint_name, $strelka_job_id, $strelka_vcf, $dirs, $opt) = @_;

    say "\n### SCHEDULING MERGE SOMATIC VCFS ###";

    my @job_ids;
    my $qsub = qsubJava($opt, "SOMVARMERGE");
    #    my $input_vcf;
    my $output_vcf = $strelka_vcf;
    my $job_id;

    #    if ($opt->{SOMVAR_TARGETS}) {
    #        $input_vcf = $strelka_vcf;
    #        $output_vcf = catfile($dirs->{out}, "${joint_name}_filtered_somatics.vcf");
    #
    #        $job_id = fromTemplate(
    #            "SomaticFiltering",
    #            undef,
    #            0,
    #            $qsub,
    #            [$strelka_job_id],
    #            $dirs,
    #            $opt,
    #            input_vcf => $input_vcf,
    #            output_vcf => $output_vcf,
    #        );
    #        push @job_ids, $job_id;
    #    }

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

    my $melted_vcf = catfile($dirs->{out}, "${joint_name}_melted_without_pon.vcf");
    $job_id = fromTemplate(
        "SomaticMelting",
        undef,
        0,
        $qsub,
        [$job_id],
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        joint_name => $joint_name,
        pre_annotate_vcf => $pre_annotate_vcf,
        input_vcf => $output_vcf,
        output_vcf => $melted_vcf,
    );
    push @job_ids, $job_id;

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
        input_vcf => $melted_vcf,
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
        joint_name => $joint_name,
        ref_bam_path => $ref_bam_path,
        tumor_bam_path => $tumor_bam_path,
        final_vcf => $final_vcf,
    );

    return ($job_id, $final_vcf);
}

1;
