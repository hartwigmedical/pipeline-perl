package HMF::Pipeline::SomaticVariants;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;
use Carp;

use HMF::Pipeline::Config qw(createDirs addSubDir);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Sge qw(qsubTemplate qsubJava);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING SOMATIC VARIANT CALLERS ###";

    my @pileup_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if ($opt->{SOMVAR_VARSCAN} eq "yes") {
            say "Creating pileup for: $sample";
            my $pileup_job = runPileup($sample, $opt);
            push @pileup_jobs, $pileup_job;
        }
    }
    $opt->{RUNNING_JOBS}->{'pileup'} = \@pileup_jobs;

    my $metadata = HMF::Pipeline::Metadata::parse($opt);
    my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
    my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";

    $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};
    $opt->{BAM_FILES}->{$tumor_sample} or die "metadata tumor_sample $tumor_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};

    my $somatic_name = "${ref_sample}_${tumor_sample}";
    my $out_dir = catfile($opt->{OUTPUT_DIR}, "somaticVariants", $somatic_name);
    my $dirs = createDirs($out_dir);

    my @running_jobs;
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$tumor_sample}} if @{$opt->{RUNNING_JOBS}->{$tumor_sample}};
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$ref_sample}} if @{$opt->{RUNNING_JOBS}->{$ref_sample}};
    my $tumor_sample_bam = catfile($opt->{OUTPUT_DIR}, $tumor_sample, "mapping", $opt->{BAM_FILES}->{$tumor_sample});
    my $ref_sample_bam = catfile($opt->{OUTPUT_DIR}, $ref_sample, "mapping", $opt->{BAM_FILES}->{$ref_sample});

    say "\n$somatic_name \t $ref_sample_bam \t $tumor_sample_bam";

    my $done_file = catfile($dirs->{log}, "${somatic_name}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my @somvar_jobs;
    if ($opt->{SOMVAR_STRELKA} eq "yes") {
        say "\n### SCHEDULING STRELKA ###";
        my $strelka_job = runStrelka($tumor_sample, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, $dirs, $opt);
        push @somvar_jobs, $strelka_job if $strelka_job;
    }

    if ($opt->{SOMVAR_VARSCAN} eq "yes") {
        say "\n### SCHEDULING VARSCAN ###";
        my $varscan_job = runVarscan($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, $dirs, $opt);
        push @somvar_jobs, $varscan_job if $varscan_job;
    }

    if ($opt->{SOMVAR_FREEBAYES} eq "yes") {
        say "\n### SCHEDULING FREEBAYES ###";
        my $freebayes_job = runFreeBayes($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, $dirs, $opt);
        push @somvar_jobs, $freebayes_job if $freebayes_job;
    }

    if ($opt->{SOMVAR_MUTECT} eq "yes") {
        say "\n### SCHEDULING MUTECT ###";
        my $mutect_job = runMutect($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, $dirs, $opt);
        push @somvar_jobs, $mutect_job if $mutect_job;
    }

    my $job_id = mergeSomatics($tumor_sample, $somatic_name, \@somvar_jobs, $dirs, $opt);
    $opt->{RUNNING_JOBS}->{somvar} = [$job_id];
    return;
}

sub mergeSomatics {
    my ($tumor_sample, $somatic_name, $somvar_jobs, $dirs, $opt) = @_;

    say "\n### SCHEDULING MERGE SOMATIC VCFS ###";

    my @inputs;
    push @inputs, "-V:strelka strelka/passed.somatic.merged.vcf" if $opt->{SOMVAR_STRELKA} eq "yes";
    push @inputs, "-V:varscan varscan/${somatic_name}.merged.Somatic.hc.vcf" if $opt->{SOMVAR_VARSCAN} eq "yes";
    push @inputs, "-V:freebayes freebayes/${somatic_name}_somatic_filtered.vcf" if $opt->{SOMVAR_FREEBAYES} eq "yes";
    push @inputs, "-V:mutect mutect/${somatic_name}_mutect_passed.vcf" if $opt->{SOMVAR_MUTECT} eq "yes";

    my $in_vcf;
    my $out_vcf = catfile($dirs->{out}, "${somatic_name}_merged_somatics.vcf");

    my $hold_jid;
    my $job_id = "SomaticMerge_${tumor_sample}_" . getId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "SomaticMerging.sh.tt", $bash_file,
        inputs => \@inputs,
        out_vcf => $out_vcf,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "SOMVARMERGE");
    if (@{$somvar_jobs}) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @{$somvar_jobs}) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }

    if ($opt->{SOMVAR_TARGETS}) {
        $in_vcf = $out_vcf;
        $out_vcf = catfile($dirs->{out}, "${somatic_name}_filtered_merged_somatics.vcf");

        $hold_jid = $job_id;
        $job_id = "SomaticFiltering_${tumor_sample}_" . getId();
        $bash_file = catfile($dirs->{job}, "${job_id}.sh");

        writeFromTemplate(
            "SomaticFiltering.sh.tt", $bash_file,
            in_vcf => $in_vcf,
            out_vcf => $out_vcf,
            dirs => $dirs,
            opt => $opt,
        );

        $qsub = qsubJava($opt, "SOMVARMERGE");
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";
    }

    my $pre_annotate_vcf = $out_vcf;

    if ($opt->{SOMVAR_ANNOTATE} eq "yes") {
        (my $basename = $out_vcf) =~ s/\.vcf$//;
        $out_vcf = "${basename}_annotated.vcf";
        $hold_jid = $job_id;
        $job_id = "SomaticAnnotate_${tumor_sample}_" . getId();
        $bash_file = catfile($dirs->{job}, "${job_id}.sh");

        writeFromTemplate(
            "SomaticAnnotation.sh.tt", $bash_file,
            basename => $basename,
            finalvcf => $out_vcf,
            dirs => $dirs,
            opt => $opt,
        );

        $qsub = qsubJava($opt, "SOMVARMERGE");
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";
    }

    $in_vcf = $out_vcf;
    $out_vcf =~ s/\.vcf$/_melted.vcf/;
    $hold_jid = $job_id;
    $job_id = "SomaticMelt_${tumor_sample}_" . getId();
    $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "SomaticMelting.sh.tt", $bash_file,
        pre_annotate_vcf => $pre_annotate_vcf,
        in_vcf => $in_vcf,
        out_vcf => $out_vcf,
        tumor_sample => $tumor_sample,
        somatic_name => $somatic_name,
        dirs => $dirs,
        opt => $opt,
    );

    $qsub = qsubJava($opt, "SOMVARMERGE");
    system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";

    HMF::Pipeline::Metadata::linkArtefact($out_vcf, "somatic_vcf", $opt);
    HMF::Pipeline::Metadata::linkArtefact("${out_vcf}.idx", "somatic_vcf_index", $opt);

    return $job_id;
}

sub runStrelka {
    my ($tumor_sample, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};

    $dirs->{strelka}->{out} = catfile($dirs->{out}, "strelka");

    my $job_id = "STR_${tumor_sample}_" . getId();
    my $done_file = catfile($dirs->{log}, "strelka.done");
    if (-f $done_file) {
        say "WARNING: $done_file, skipping $job_id";
        return;
    }

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "Strelka.sh.tt", $bash_file,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "STRELKA");
    if (@running_jobs) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }

    return $job_id;
}

sub runPileup {
    my ($sample, $opt) = @_;

    my $bam = $opt->{BAM_FILES}->{$sample};
    (my $pileup = $bam) =~ s/\.bam/\.pileup/;

    my $job_id = "PileUp_${sample}_" . getId();
    my $done_file = catfile($opt->{OUTPUT_DIR}, $sample, "logs", "Pileup_${sample}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return $job_id;
    }

    my $log_dir = catfile($opt->{OUTPUT_DIR}, $sample, "logs");
    my $bash_file = catfile($opt->{OUTPUT_DIR}, $sample, "jobs", "${job_id}.sh");

    writeFromTemplate(
        "PileUp.sh.tt", $bash_file,
        sample => $sample,
        bam => $bam,
        pileup => $pileup,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "PILEUP");
    if (@{$opt->{RUNNING_JOBS}->{$sample}}) {
        system "$qsub -o $log_dir/Pileup_$sample.out -e $log_dir/Pileup_$sample.err -N $job_id -hold_jid " . join(",", @{$opt->{RUNNING_JOBS}->{$sample}}) . " $bash_file";
    } else {
        system "$qsub -o $log_dir/Pileup_$sample.out -e $log_dir/Pileup_$sample.err -N $job_id $bash_file";
    }
    return $job_id;
}

sub runVarscan {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{'pileup'}};

    $dirs->{varscan}->{out} = addSubDir($dirs, "varscan");

    my $job_id = "VS_${tumor_sample}_" . getId();
    my $done_file = catfile($dirs->{log}, "varscan.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    (my $tumor_sample_pileup = $tumor_sample_bam) =~ s/\.bam$/\.pileup\.gz/;
    (my $ref_sample_pileup = $ref_sample_bam) =~ s/\.bam$/\.pileup\.gz/;

    my @chrs = @{getChromosomes($opt)};
    my @varscan_jobs;
    foreach my $chr (@chrs) {
        my $job_id = "VS_${tumor_sample}_${chr}_" . getId();
        my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        writeFromTemplate(
            "Varscan.sh.tt", $bash_file,
            chr => $chr,
            output_name => $output_name,
            ref_sample_pileup => $ref_sample_pileup,
            tumor_sample_pileup => $tumor_sample_pileup,
            dirs => $dirs,
            opt => $opt,
        );

        my $qsub = qsubJava($opt, "VARSCAN");
        if (@running_jobs) {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
        }

        push @varscan_jobs, $job_id;
    }

    my @snp_vcfs = map { "${somatic_name}_${_}.snp.vcf" } @chrs;
    my @indel_vcfs = map { "${somatic_name}_${_}.indel.vcf" } @chrs;

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "VarscanPS.sh.tt", $bash_file,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        snp_vcfs => \@snp_vcfs,
        indel_vcfs => \@indel_vcfs,
        ref_sample_pileup => $ref_sample_pileup,
        tumor_sample_pileup => $tumor_sample_pileup,
        somatic_name => $somatic_name,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "VARSCAN");
    if (@varscan_jobs) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @varscan_jobs) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

sub runFreeBayes {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};

    $dirs->{freebayes}->{out} = addSubDir($dirs, "freebayes");
    $dirs->{freebayes}->{tmp} = addSubDir($dirs->{freebayes}, "tmp");

    my $job_id = "FB_${tumor_sample}_" . getId();
    my $done_file = catfile($dirs->{log}, "freebayes.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    my @chrs = @{getChromosomes($opt)};
    my @freebayes_jobs;
    foreach my $chr (@chrs) {
        my $job_id = "FB_${tumor_sample}_${chr}_" . getId();
        my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        writeFromTemplate(
            "Freebayes.sh.tt", $bash_file,
            chr => $chr,
            tumor_sample_bam => $tumor_sample_bam,
            ref_sample_bam => $ref_sample_bam,
            output_name => $output_name,
            dirs => $dirs,
            opt => $opt,
        );

        my $qsub = qsubJava($opt, "FREEBAYES");
        if (@running_jobs) {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
        }

        push @freebayes_jobs, $job_id;
    }

    my @snp_vcfs = map { "${somatic_name}_${_}.vcf" } @chrs;

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "FreebayesPostProcess.sh.tt", $bash_file,
        snp_vcfs => \@snp_vcfs,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        somatic_name => $somatic_name,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "FREEBAYES");
    if (@freebayes_jobs) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @freebayes_jobs) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

sub runMutect {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};

    $dirs->{mutect}->{out} = addSubDir($dirs, "mutect");
    $dirs->{mutect}->{tmp} = addSubDir($dirs->{mutect}, "tmp");

    my $job_id = "MUT_${tumor_sample}_" . getId();
    my $done_file = catfile($dirs->{log}, "mutect.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    my @chrs = @{getChromosomes($opt)};
    my @mutect_jobs;
    foreach my $chr (@chrs) {
        my $job_id = "MUT_${tumor_sample}_${chr}_" . getId();
        my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        writeFromTemplate(
            "Mutect.sh.tt", $bash_file,
            chr => $chr,
            tumor_sample_bam => $tumor_sample_bam,
            ref_sample_bam => $ref_sample_bam,
            output_name => $output_name,
            dirs => $dirs,
            opt => $opt,
        );

        my $qsub = qsubJava($opt, "MUTECT");
        if (@running_jobs) {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
        }
        push @mutect_jobs, $job_id;
    }

    my @vcfs = map { "${somatic_name}_${_}_mutect.vcf" } @chrs;

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "MutectCF.sh.tt", $bash_file,
        vcfs => \@vcfs,
        ref_sample_bam => $ref_sample_bam,
        tumor_sample_bam => $tumor_sample_bam,
        somatic_name => $somatic_name,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubJava($opt, "MUTECT");
    if (@mutect_jobs) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @mutect_jobs) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }

    return $job_id;
}

sub getChromosomes {
    my ($opt) = @_;

    (my $dict_file = $opt->{GENOME}) =~ s/\.fasta$/.dict/;
    my @chrs;
    open my $fh, "<", $dict_file or confess "could not open $dict_file: $!";
    while (<$fh>) {
        chomp;
        push @chrs, $1 if /SN:(\w+)\s*LN:(\d+)/;
    }
    close $fh;

    return \@chrs;
}

1;
