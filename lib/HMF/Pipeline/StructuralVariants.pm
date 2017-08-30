package HMF::Pipeline::StructuralVariants;

use FindBin::libs;
use discipline;

use File::Spec::Functions;
use Sort::Key::Natural qw(mkkey_natural);

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs sampleBamsAndJobs sampleControlBamsAndJobs getChromosomes);
use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile);
use HMF::Pipeline::Job::Vcf;
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Metadata qw(linkExtraArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SV CALLING ###";

    $opt->{RUNNING_JOBS}->{'sv'} = [];
    if ($opt->{SV_DELLY} eq "yes") {
        my $delly_jobs = runDelly($opt);
        push @{$opt->{RUNNING_JOBS}->{'sv'}}, @{$delly_jobs};
    }
    if ($opt->{SV_MANTA} eq "yes") {
        my $manta_jobs = runManta($opt);
        push @{$opt->{RUNNING_JOBS}->{'sv'}}, @{$manta_jobs};
    }
    return;
}

sub runDelly {
    my ($opt) = @_;

    say "\n### SCHEDULING DELLY ###";

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "delly"));
    my ($sample_bams, $running_jobs) = sampleBamsAndJobs($opt);

    my @sv_types = split /\t/, $opt->{DELLY_SVTYPE};
    my $delly_job_ids = [];
    foreach my $sv_type (@sv_types) {
        my $concat_vcf = catfile($dirs->{tmp}, "$opt->{RUN_NAME}_${sv_type}.vcf");
        my $output_vcf = catfile($dirs->{out}, "$opt->{RUN_NAME}_${sv_type}.vcf");
        my $done_file = checkReportedDoneFile("Delly", $sv_type, $dirs, $opt) or next;

        my ($chunk_job_ids, $vcf_files) = runDellyType($sample_bams, $sv_type, $running_jobs, $dirs, $opt);
        my $concat_job_id = HMF::Pipeline::Job::Vcf::concat($vcf_files, $concat_vcf, $sv_type, "DELLY_MERGE", $chunk_job_ids, $dirs, $opt);
        my $sort_job_id = HMF::Pipeline::Job::Vcf::sorted($concat_vcf, $output_vcf, $sv_type, "DELLY_MERGE", [$concat_job_id], $dirs, $opt);
        my $job_id = HMF::Pipeline::Job::markDone($done_file, [ @{$chunk_job_ids}, $concat_job_id, $sort_job_id ], $dirs, $opt);
        push @{$delly_job_ids}, $job_id;
    }
    return $delly_job_ids;
}

sub runDellyType {
    my ($sample_bams, $sv_type, $running_jobs, $dirs, $opt) = @_;

    # DEL, DUP, INS, INV are only ever intra-chromosomal in Delly
    if ($sv_type eq "TRA") {
        return runDellyInterchromosomal($sample_bams, $sv_type, $running_jobs, $dirs, $opt);
    } else {
        return runDellyIntrachromosomal($sample_bams, $sv_type, $running_jobs, $dirs, $opt);
    }
}

sub runDellyInterchromosomal {
    my ($sample_bams, $type, $running_jobs, $dirs, $opt) = @_;

    my ($job_ids, $vcf_files) = ([], []);
    my @chrs = @{getChromosomes($opt)};
    foreach my $chr1 (@chrs) {
        foreach my $chr2 (@chrs) {
            next if mkkey_natural($chr1) le mkkey_natural($chr2);
            my ($job_id, $output_vcf) = runDellyJob(
                "${type}_${chr1}_${chr2}",
                # include these
                {$chr1 => 1, $chr2 => 1},
                $job_ids,
                $vcf_files,
                $sample_bams,
                $type,
                $running_jobs,
                $dirs,
                $opt,
            );
        }
    }
    return ($job_ids, $vcf_files);
}

sub runDellyIntrachromosomal {
    my ($sample_bams, $type, $running_jobs, $dirs, $opt) = @_;

    my ($job_ids, $vcf_files) = ([], []);
    foreach my $chr (@{getChromosomes($opt)}) {
        my ($job_id, $output_vcf) = runDellyJob(
            "${type}_${chr}",
            # include these
            {$chr => 1},
            $job_ids,
            $vcf_files,
            $sample_bams,
            $type,
            $running_jobs,
            $dirs,
            $opt,
        );
    }
    return ($job_ids, $vcf_files);
}

sub runDellyJob {
    my ($step, $include_chrs, $job_ids, $vcf_files, $sample_bams, $type, $running_jobs, $dirs, $opt) = @_;

    my $exclude_file = catfile($dirs->{tmp}, "${step}_exclude.txt");
    open my $fh, ">", $exclude_file or die "could not open $exclude_file: $!";
    foreach my $exclude_chr (@{getChromosomes($opt)}) {
        say $fh $exclude_chr unless defined $include_chrs->{$exclude_chr};
    }
    close $fh;

    my $output_vcf = catfile($dirs->{tmp}, "${step}.vcf");
    my $job_id = fromTemplate(
        "Delly",
        $step,
        0,
        qsubTemplate($opt, "DELLY"),
        $running_jobs,
        $dirs,
        $opt,
        step => $step,
        sample_bams => $sample_bams,
        type => $type,
        exclude_file => $exclude_file,
        output_vcf => $output_vcf,
    );

    push @{$job_ids}, $job_id;
    push @{$vcf_files}, $output_vcf;
    return;
}

sub runManta {
    my ($opt) = @_;

    say "\n### SCHEDULING MANTA ###";

    if ($opt->{SV_MODE} eq "sample") {
        my @manta_jobs;
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);
            say "\n$sample \t $sample_bam";
            my $job_id = runMantaJob($sample_bam, undef, $sample, $running_jobs, $opt);
            push @manta_jobs, $job_id;
        }
        return \@manta_jobs;
    } else {
        my @manta_jobs;
        my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
        say "\n$joint_name \t $ref_sample_bam \t $tumor_sample_bam";
        my $job_id = runMantaJob($tumor_sample_bam, $ref_sample_bam, $joint_name, $running_jobs, $opt);
        push @manta_jobs, $job_id;
        $job_id = runBreakpointInspector($tumor_sample, $tumor_sample_bam, $ref_sample, $ref_sample_bam, $joint_name, $job_id, $opt);
        push @manta_jobs, $job_id;
        return \@manta_jobs;
    }
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

    linkExtraArtefact($opt->{BPI_VCF_FILE}, $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${control}_sliced.bam"), $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${sample}_sliced.bam"), $opt);
    linkExtraArtefact(catfile($dirs->{out}, "${joint_name}_bpi_stats.tsv"), $opt);

    return $job_id;
}

1;