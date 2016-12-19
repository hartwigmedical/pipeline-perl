package HMF::Pipeline::StructuralVariants;

use FindBin::libs;
use discipline;

use File::Spec::Functions;
use Sort::Key::Natural qw(mkkey_natural);

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs sampleBamsAndJobs sampleControlBamsAndJobs getChromosomes);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Job::Vcf;
use HMF::Pipeline::Sge qw(qsubTemplate);

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

    my %sv_actions;
    my @sv_types = split /\t/, $opt->{DELLY_SVTYPE};
    my @sv_split = map { $_ eq "yes" } split /\t/, $opt->{DELLY_SPLIT};
    @sv_actions{@sv_types} = @sv_split;

    my $delly_job_ids = [];
    while (my ($type, $split) = each %sv_actions) {
        my $concat_vcf = catfile($dirs->{tmp}, "$opt->{RUN_NAME}_${type}.vcf");
        my $output_vcf = catfile($dirs->{out}, "$opt->{RUN_NAME}_${type}.vcf");
        my $concat_job_id;

        if ($split) {
            # DEL, INF, DUP are only ever intra-chromosomal in Delly
            my ($chunk_job_ids, $vcf_files);
            if ($type eq "TRA") {
                ($chunk_job_ids, $vcf_files) = runDellyInterchromosomal($sample_bams, $type, $running_jobs, $dirs, $opt);
            } else {
                ($chunk_job_ids, $vcf_files) = runDellyIntrachromosomal($sample_bams, $type, $running_jobs, $dirs, $opt);
            }
            $concat_job_id = HMF::Pipeline::Job::Vcf::concat($vcf_files, $concat_vcf, $type, "DELLY_MERGE", $chunk_job_ids, $dirs, $opt);
        } else {
            $concat_job_id = runDellyJob($type, $sample_bams, $type, undef, $concat_vcf, $running_jobs, $dirs, $opt);
        }
        my $sort_job_id = HMF::Pipeline::Job::Vcf::sorted($concat_vcf, $output_vcf, $type, "DELLY_MERGE", [$concat_job_id], $dirs, $opt);
        push @{$delly_job_ids}, $sort_job_id if $sort_job_id;
    }
    return $delly_job_ids;
}

sub runDellyInterchromosomal {
    my ($sample_bams, $type, $running_jobs, $dirs, $opt) = @_;

    my (@job_ids, @vcf_files);
    my @chrs = @{getChromosomes($opt)};
    foreach my $chr1 (@chrs) {
        foreach my $chr2 (@chrs) {
            next unless mkkey_natural($chr2) gt mkkey_natural($chr1);
            my $step = "${type}_${chr1}_${chr2}";

            my $exclude_file = catfile($dirs->{tmp}, "${step}_exclude.txt");
            open my $fh, ">", $exclude_file;
            foreach my $exclude_chr (@chrs) {
                say $fh $exclude_chr unless $exclude_chr eq $chr1 or $exclude_chr eq $chr2;
            }
            close $fh;

            my $output_vcf = catfile($dirs->{tmp}, "${step}.vcf");
            my $job_id = runDellyJob($step, $sample_bams, $type, $exclude_file, $output_vcf, $running_jobs, $dirs, $opt);
            push @job_ids, $job_id;
            push @vcf_files, $output_vcf;
        }
    }
    return (\@job_ids, \@vcf_files);
}

sub runDellyIntrachromosomal {
    my ($sample_bams, $type, $running_jobs, $dirs, $opt) = @_;

    my (@job_ids, @vcf_files);
    my @chrs = @{getChromosomes($opt)};
    foreach my $chr (@chrs) {
        my $step = "${type}_${chr}";

        my $exclude_file = catfile($dirs->{tmp}, "${step}_exclude.txt");
        open my $fh, ">", $exclude_file;
        foreach my $exclude_chr (@chrs) {
            say $fh $exclude_chr unless $exclude_chr eq $chr;
        }
        close $fh;

        my $output_vcf = catfile($dirs->{tmp}, "${step}.vcf");
        my $job_id = runDellyJob($step, $sample_bams, $type, $exclude_file, $output_vcf, $running_jobs, $dirs, $opt);
        push @job_ids, $job_id;
        push @vcf_files, $output_vcf;
    }
    return (\@job_ids, \@vcf_files);
}

sub runDellyJob {
    my ($step, $sample_bams, $type, $exclude_file, $output_vcf, $running_jobs, $dirs, $opt) = @_;

    my $job_id = fromTemplate(
        "Delly",
        $step,
        qsubTemplate($opt, "DELLY"),
        $running_jobs,
        $dirs,
        $opt,
        step => $step,
        sample_bams => [ values %{$sample_bams} ],
        type => $type,
        exclude_file => $exclude_file,
        output_vcf => $output_vcf,
    );

    return $job_id;
}

sub runManta {
    my ($opt) = @_;

    say "\n### SCHEDULING MANTA ###";

    if ($opt->{SV_MODE} eq "sample") {
        my @manta_jobs;
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my ($sample_bam, $running_jobs) = sampleBamAndJobs($sample, $opt);
            say "\n$sample \t $sample_bam";
            my $job_id = runMantaJob($sample, $sample_bam, undef, undef, $sample, $running_jobs, $opt);
            next if not $job_id;
            push @manta_jobs, $job_id;
        }
        return \@manta_jobs;
    } else {
        my ($ref_sample, $tumor_sample, $ref_sample_bam, $tumor_sample_bam, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
        say "\n$joint_name \t $ref_sample_bam \t $tumor_sample_bam";
        my $job_id = runMantaJob($tumor_sample, $tumor_sample_bam, $ref_sample, $ref_sample_bam, $joint_name, $running_jobs, $opt);
        return [$job_id];
    }
}

sub runMantaJob {
    my ($sample, $sample_bam, $control, $control_bam, $joint_name, $running_jobs, $opt) = @_;

    my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, "structuralVariants", "manta", $joint_name));
    my $job_id = fromTemplate(
        "Manta",
        $joint_name,
        qsubTemplate($opt, "MANTA"),
        $running_jobs,
        $dirs,
        $opt,
        sample => $sample,
        sample_bam => $sample_bam,
        control => $control,
        control_bam => $control_bam,
        joint_name => $joint_name,
    );

    return $job_id;
}

1;
