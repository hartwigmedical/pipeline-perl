package HMF::Pipeline::SomaticVariants;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs addSubDir getChromosomes sampleControlBamsAndJobs);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Metadata;
use HMF::Pipeline::Sge qw(qsubTemplate qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SOMATIC VARIANT CALLERS ###";

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "somaticVariants", $joint_name));
    say "\n$joint_name \t $ref_bam_path \t $tumor_bam_path";

    # must be before returning if already .done
    my $final_vcf = catfile($dirs->{out}, "${joint_name}_melted.vcf");
    HMF::Pipeline::Metadata::linkArtefact($final_vcf, "somatic_vcf", $opt);
    HMF::Pipeline::Metadata::linkArtefact("${final_vcf}.idx", "somatic_vcf_index", $opt);

    my $done_file = checkReportedDoneFile($joint_name, undef, $dirs, $opt) or return;

    my (@somvar_jobs, %somvar_vcfs);
    if ($opt->{SOMVAR_FREEBAYES} eq "yes") {
        my ($job_id, $vcf) = runFreebayes($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt);
        push @somvar_jobs, $job_id;
        $somvar_vcfs{freebayes} = $vcf;
    }
    if ($opt->{SOMVAR_MUTECT} eq "yes") {
        my ($job_id, $vcf) = runMutect($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt);
        push @somvar_jobs, $job_id;
        $somvar_vcfs{mutect} = $vcf;
    }
    if ($opt->{SOMVAR_STRELKA} eq "yes") {
        my ($job_id, $vcf) = runStrelka($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt);
        push @somvar_jobs, $job_id;
        $somvar_vcfs{strelka} = $vcf;
    }
    if ($opt->{SOMVAR_VARSCAN} eq "yes") {
        my ($job_id, $vcf) = runVarscan($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt);
        push @somvar_jobs, $job_id;
        $somvar_vcfs{varscan} = $vcf;
    }

    my $merge_job_id = mergeSomatics($tumor_sample, $joint_name, \@somvar_jobs, \%somvar_vcfs, $final_vcf, $dirs, $opt);
    my $job_id = markDone($done_file, [ @somvar_jobs, $merge_job_id ], $dirs, $opt);
    $opt->{RUNNING_JOBS}->{somvar} = [$job_id];

    return;
}

sub mergeSomatics {
    my ($tumor_sample, $joint_name, $somvar_jobs, $somvar_vcfs, $final_vcf, $dirs, $opt) = @_;

    say "\n### SCHEDULING MERGE SOMATIC VCFS ###";

    my $in_vcf;
    my $out_vcf = catfile($dirs->{out}, "${joint_name}_merged_somatics.vcf");
    my $job_id = fromTemplate(
        "SomaticMerging",
        undef,
        0,
        qsubJava($opt, "SOMVARMERGE"),
        $somvar_jobs,
        $dirs,
        $opt,
        in_vcfs => $somvar_vcfs,
        out_vcf => $out_vcf,
    );

    if ($opt->{SOMVAR_TARGETS}) {
        $in_vcf = $out_vcf;
        $out_vcf = catfile($dirs->{out}, "${joint_name}_filtered_merged_somatics.vcf");

        $job_id = fromTemplate(
            "SomaticFiltering",
            undef,
            0,
            qsubJava($opt, "SOMVARMERGE"),
            [$job_id],
            $dirs,
            $opt,
            in_vcf => $in_vcf,
            out_vcf => $out_vcf,
        );
    }

    my $pre_annotate_vcf = $out_vcf;
    if ($opt->{SOMVAR_ANNOTATE} eq "yes") {
        (my $basename = $out_vcf) =~ s/\.vcf$//;
        $out_vcf = "${basename}_annotated.vcf";

        $job_id = fromTemplate(
            "SomaticAnnotation",
            undef,
            0,
            qsubJava($opt, "SOMVARMERGE"),
            [$job_id],
            $dirs,
            $opt,
            basename => $basename,
            final_vcf => $out_vcf,
        );
    }

    $job_id = fromTemplate(
        "SomaticMelting",
        undef,
        0,
        qsubJava($opt, "SOMVARMERGE"),
        [$job_id],
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        joint_name => $joint_name,
        pre_annotate_vcf => $pre_annotate_vcf,
        in_vcf => $out_vcf,
        out_vcf => $final_vcf,
    );

    return $job_id;
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

sub runVarscan {
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING VARSCAN ###";

    $dirs->{varscan}->{out} = addSubDir($dirs, "varscan");
    my $final_vcf = catfile($dirs->{varscan}->{out}, "${joint_name}.merged.Somatic.hc.vcf");

    my $done_file = checkReportedDoneFile("Varscan", undef, $dirs, $opt) or return (undef, $final_vcf);
    my (@chr_jobs, @chr_snp_vcfs, @chr_indel_vcfs);
    foreach my $chr (@{getChromosomes($opt)}) {
        my $snp_vcf = "${joint_name}_${chr}.snp.vcf";
        my $indel_vcf = "${joint_name}_${chr}.indel.vcf";

        my $job_id = fromTemplate(
            "Varscan",
            $chr,
            0,
            qsubJava($opt, "VARSCAN"),
            [ @{$running_jobs}, @{$opt->{RUNNING_JOBS}->{pileup}} ],
            $dirs,
            $opt,
            chr => $chr,
            ref_pileup => $opt->{PILEUP_FILES}->{$ref_sample},
            tumor_pileup => $opt->{PILEUP_FILES}->{$tumor_sample},
            snp_vcf => $snp_vcf,
            indel_vcf => $indel_vcf,
        );
        push @chr_jobs, $job_id;
        push @chr_snp_vcfs, $snp_vcf;
        push @chr_indel_vcfs, $indel_vcf;
    }

    my $post_job_id = fromTemplate(
        "VarscanPostProcess",
        undef,
        0,
        qsubJava($opt, "VARSCAN"),
        \@chr_jobs,
        $dirs,
        $opt,
        joint_name => $joint_name,
        snp_vcfs => \@chr_snp_vcfs,
        indel_vcfs => \@chr_indel_vcfs,
        final_vcf => $final_vcf,
    );
    my $job_id = markDone($done_file, [ @chr_jobs, $post_job_id ], $dirs, $opt);
    return ($job_id, $final_vcf);
}

sub runFreebayes {
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING FREEBAYES ###";

    $dirs->{freebayes}->{out} = addSubDir($dirs, "freebayes");
    $dirs->{freebayes}->{tmp} = addSubDir($dirs->{freebayes}, "tmp");
    my $final_vcf = catfile($dirs->{freebayes}->{out}, "${joint_name}_somatic_filtered.vcf");

    my $done_file = checkReportedDoneFile("Freebayes", undef, $dirs, $opt) or return (undef, $final_vcf);
    my (@chr_jobs, @chr_vcfs);
    foreach my $chr (@{getChromosomes($opt)}) {
        my $output_vcf = catfile($dirs->{freebayes}->{out}, "${joint_name}_${chr}.vcf");
        my $job_id = fromTemplate(
            "Freebayes",
            $chr,
            0,
            qsubJava($opt, "FREEBAYES"),
            $running_jobs,
            $dirs,
            $opt,
            chr => $chr,
            ref_bam_path => $ref_bam_path,
            tumor_bam_path => $tumor_bam_path,
            output_vcf => $output_vcf,
        );
        push @chr_jobs, $job_id;
        push @chr_vcfs, $output_vcf;
    }

    my $post_job_id = fromTemplate(
        "FreebayesPostProcess",
        undef,
        0,
        qsubJava($opt, "FREEBAYES"),
        \@chr_jobs,
        $dirs,
        $opt,
        joint_name => $joint_name,
        input_vcfs => \@chr_vcfs,
        final_vcf => $final_vcf,
    );
    my $job_id = markDone($done_file, [ @chr_jobs, $post_job_id ], $dirs, $opt);
    return ($job_id, $final_vcf);
}

sub runMutect {
    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING MUTECT ###";

    $dirs->{mutect}->{out} = addSubDir($dirs, "mutect");
    $dirs->{mutect}->{tmp} = addSubDir($dirs->{mutect}, "tmp");
    my $final_vcf = catfile($dirs->{mutect}->{out}, "${joint_name}_mutect_passed.vcf");

    my $done_file = checkReportedDoneFile("Mutect", undef, $dirs, $opt) or return (undef, $final_vcf);
    my (@chr_jobs, @chr_vcfs);
    foreach my $chr (@{getChromosomes($opt)}) {
        my $output_vcf = catfile($dirs->{mutect}->{tmp}, "${joint_name}_${chr}.vcf");
        my $job_id = fromTemplate(
            "Mutect",
            $chr,
            0,
            qsubJava($opt, "MUTECT"),
            $running_jobs,
            $dirs,
            $opt,
            chr => $chr,
            ref_bam_path => $ref_bam_path,
            tumor_bam_path => $tumor_bam_path,
            output_vcf => $output_vcf,
        );
        push @chr_jobs, $job_id;
        push @chr_vcfs, $output_vcf;
    }

    my $post_job_id = fromTemplate(
        "MutectPostProcess",
        undef,
        0,
        qsubJava($opt, "MUTECT"),
        \@chr_jobs,
        $dirs,
        $opt,
        joint_name => $joint_name,
        input_vcfs => \@chr_vcfs,
        final_vcf => $final_vcf,
    );
    my $job_id = markDone($done_file, [ @chr_jobs, $post_job_id ], $dirs, $opt);
    return ($job_id, $final_vcf);
}

1;
