package illumina_somaticVariants;

use 5.16.0;
use strict;
use warnings;

use File::Path qw(make_path);
use File::Basename;
use File::Spec::Functions;
use Carp;

use FindBin;
use lib "$FindBin::Bin";

use illumina_sge qw(getJobId qsubTemplate qsubJava);
use illumina_template qw(from_template);
use illumina_metadataParser qw(metadataParse);


sub runSomaticVariantCallers {
    my ($opt) = @_;

    my @pileupJobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        if ($opt->{SOMVAR_VARSCAN} eq "yes") {
            say "Creating pileup for: $sample";
            my $pileup_job = runPileup($sample, $opt);
            push @pileupJobs, $pileup_job;
        }
    }

    $opt->{RUNNING_JOBS}->{'pileup'} = \@pileupJobs;

    my $metadata = metadataParse($opt->{OUTPUT_DIR});

    my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
    my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";

    $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list";
    $opt->{BAM_FILES}->{$tumor_sample} or die "metadata tumor_sample $tumor_sample not in BAM file list";

    my $somatic_name = "${ref_sample}_${tumor_sample}";
    my $out_dir = catfile($opt->{OUTPUT_DIR}, "somaticVariants", $somatic_name);
    my %somatic_dirs = (
        out => $out_dir,
        tmp => catfile($out_dir, "tmp"),
        log => catfile($out_dir, "logs"),
        job => catfile($out_dir, "jobs"),
    );

    make_path(values %somatic_dirs, { error => \my $errors });
    my $messages = join ", ", map { join ": ", each $_ } @{$errors};
    die "Couldn't create somatic output directories: $messages" if $messages;

    my @running_jobs;
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$tumor_sample}} if @{$opt->{RUNNING_JOBS}->{$tumor_sample}};
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$ref_sample}} if @{$opt->{RUNNING_JOBS}->{$ref_sample}};
    my $tumor_sample_bam = catfile($opt->{OUTPUT_DIR}, $tumor_sample, "mapping", $opt->{BAM_FILES}->{$tumor_sample});
    my $ref_sample_bam = catfile($opt->{OUTPUT_DIR}, $ref_sample, "mapping", $opt->{BAM_FILES}->{$ref_sample});

    say "\n$somatic_name \t $ref_sample_bam \t $tumor_sample_bam";

    my $done_file = catfile($somatic_dirs{log}, "${somatic_name}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my @somvar_jobs;
    if ($opt->{SOMVAR_STRELKA} eq "yes") {
        say "\n###SCHEDULING STRELKA####";
        my $strelka_job = runStrelka($tumor_sample, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, \%somatic_dirs, $opt);
        push @somvar_jobs, $strelka_job if $strelka_job;
    }

    if ($opt->{SOMVAR_VARSCAN} eq "yes") {
        say "\n###SCHEDULING VARSCAN####";
        my $varscan_job = runVarscan($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, \%somatic_dirs, $opt);
        push @somvar_jobs, $varscan_job if $varscan_job;
    }

    if ($opt->{SOMVAR_FREEBAYES} eq "yes") {
        say "\n###SCHEDULING FREEBAYES####";
        my $freebayes_job = runFreeBayes($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, \%somatic_dirs, $opt);
        push @somvar_jobs, $freebayes_job if $freebayes_job;
    }

    if ($opt->{SOMVAR_MUTECT} eq "yes") {
        say "\n###SCHEDULING MUTECT####";
        my $mutect_job = runMutect($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, \@running_jobs, \%somatic_dirs, $opt);
        push @somvar_jobs, $mutect_job if $mutect_job;
    }

    my $job_id = mergeSomatics($tumor_sample, $somatic_name, \@somvar_jobs, \%somatic_dirs, $opt);
    $opt->{RUNNING_JOBS}->{somvar} = [$job_id];
    return;
}

sub mergeSomatics {
    my ($tumor_sample, $somatic_name, $somvar_jobs, $somatic_dirs, $opt) = @_;
    my $runName = basename($opt->{OUTPUT_DIR});

    say "\n###SCHEDULING MERGE SOMATIC VCFS####";

    my @inputs;
    push @inputs, "-V:strelka strelka/passed.somatic.merged.vcf" if $opt->{SOMVAR_STRELKA} eq "yes";
    push @inputs, "-V:varscan varscan/${somatic_name}.merged.Somatic.hc.vcf" if $opt->{SOMVAR_VARSCAN} eq "yes";
    push @inputs, "-V:freebayes freebayes/${somatic_name}_somatic_filtered.vcf" if $opt->{SOMVAR_FREEBAYES} eq "yes";
    push @inputs, "-V:mutect mutect/${somatic_name}_mutect_passed.vcf" if $opt->{SOMVAR_MUTECT} eq "yes";

    my $invcf;
    my $outvcf = catfile($somatic_dirs->{out}, "${somatic_name}_merged_somatics.vcf");

    my $hold_jid;
    my $job_id = "SomaticMerge_${tumor_sample}_" . getJobId();
    my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    from_template("SomaticMerging.sh.tt", $bash_file,
                  inputs => join(" ", @inputs),
                  outvcf => $outvcf,
                  dirs => $somatic_dirs,
                  opt => $opt,
                  runName => $runName);

    my $qsub = qsubJava($opt, "SOMVARMERGE");
    if (@{$somvar_jobs}) {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @{$somvar_jobs}) . " $bash_file";
    } else {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
    }

    if ($opt->{SOMVAR_TARGETS}) {
        $invcf = $outvcf;
        $outvcf = catfile($somatic_dirs->{out}, "${somatic_name}_filtered_merged_somatics.vcf");

        $hold_jid = $job_id;
        $job_id = "SomaticFilter_${tumor_sample}_" . getJobId();
        $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

        from_template("SomaticFiltering.sh.tt", $bash_file,
                      invcf => $invcf,
                      outvcf => $outvcf,
                      dirs => $somatic_dirs,
                      opt => $opt,
                      runName => $runName);

        $qsub = qsubJava($opt, "SOMVARMERGE");
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";
    }

    my $pre_annotate_vcf = $outvcf;

    if ($opt->{SOMVAR_ANNOTATE} eq "yes") {
        (my $basename = $outvcf) =~ s/\.vcf$//;
        $outvcf = "${basename}_annotated.vcf";
        $hold_jid = $job_id;
        $job_id = "SomaticAnnotate_${tumor_sample}_" . getJobId();
        $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

        from_template("SomaticAnnotation.sh.tt", $bash_file,
                      basename => $basename,
                      finalvcf => $outvcf,
                      dirs => $somatic_dirs,
                      opt => $opt,
                      runName => $runName);

        $qsub = qsubJava($opt, "SOMVARMERGE");
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";
    }

    $invcf = $outvcf;
    $outvcf =~ s/\.vcf$/_melted.vcf/;
    $hold_jid = $job_id;
    $job_id = "SomaticMelt_${tumor_sample}_" . getJobId();
    $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    from_template("SomaticMelting.sh.tt", $bash_file,
                  pre_annotate_vcf => $pre_annotate_vcf,
                  invcf => $invcf,
                  outvcf => $outvcf,
                  tumor_sample => $tumor_sample,
                  somatic_name => $somatic_name,
                  dirs => $somatic_dirs,
                  opt => $opt,
                  runName => $runName);

    $qsub = qsubJava($opt, "SOMVARMERGE");
    system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid $hold_jid $bash_file";

    return $job_id;
}

sub runStrelka {
    my ($tumor_sample, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $somatic_dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $runName = basename($opt{OUTPUT_DIR});

    my $done_file = catfile($somatic_dirs->{log}, "strelka.done");
    if (-f $done_file) {
        say "WARNING: $done_file, skipping";
        return;
    }

    my $job_id = "STR_".$tumor_sample."_" . getJobId();
    my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    from_template("Strelka.sh.tt", $bash_file,
                  ref_sample_bam => $ref_sample_bam,
                  tumor_sample_bam => $tumor_sample_bam,
                  dirs => $somatic_dirs,
                  runName => $runName,
                  opt => \%opt);

    my $qsub = qsubJava(\%opt, "STRELKA");
    if (@running_jobs) {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
    }

    return $job_id;
}

sub runPileup {
    my ($sample, $configuration) = @_;
    my %opt = %{$configuration};
    my $runName = basename($opt{OUTPUT_DIR});

    my $bam = $opt{BAM_FILES}->{$sample};
    (my $pileup = $bam) =~ s/\.bam/\.pileup/;
    my $jobID = "PileUp_${sample}_" . getJobId();

    my $done_file = catfile($opt{OUTPUT_DIR}, $sample, "logs", "Pileup_${sample}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return $jobID;
    }

    my $logDir = catfile($opt{OUTPUT_DIR}, $sample, "logs");
    my $bashFile = catfile($opt{OUTPUT_DIR}, $sample, "jobs", "${jobID}.sh");

    from_template("PileUp.sh.tt", "$bashFile",
                  sample => $sample,
                  bam => $bam,
                  pileup => $pileup,
                  runName => $runName,
                  opt => $configuration);

    my $qsub = qsubTemplate(\%opt, "PILEUP");
    if (@{$opt{RUNNING_JOBS}->{$sample}}) {
        system "$qsub -o $logDir/Pileup_$sample.out -e $logDir/Pileup_$sample.err -N $jobID -hold_jid " . join(",", @{$opt{RUNNING_JOBS}->{$sample}}) . " $bashFile";
    } else {
        system "$qsub -o $logDir/Pileup_$sample.out -e $logDir/Pileup_$sample.err -N $jobID $bashFile";
    }
    return $jobID;
}

sub runVarscan {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $somatic_dirs, $opt) = @_;
    my %opt = %{$opt};
    my @running_jobs = @{$running_jobs};
    push @running_jobs, @{$opt{RUNNING_JOBS}->{'pileup'}};

    my $varscan_out_dir = catfile($somatic_dirs->{out}, "varscan");
    (my $tumor_sample_pileup = $tumor_sample_bam) =~ s/\.bam$/\.pileup\.gz/;
    (my $ref_sample_pileup = $ref_sample_bam) =~ s/\.bam$/\.pileup\.gz/;

    if (!-d $varscan_out_dir) {
        make_path($varscan_out_dir) or die "Couldn't create directory $varscan_out_dir: $!";
    }

    my $done_file = catfile($somatic_dirs->{log}, "varscan.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $runName = basename($opt{OUTPUT_DIR});
    my @chrs = @{getChromosomes($opt)};
    my @varscan_jobs;

    foreach my $chr (@chrs) {
        my $job_id = "VS_${tumor_sample}_${chr}_" . getJobId();
        my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        from_template("Varscan.sh.tt", $bash_file,
                      chr => $chr,
                      output_name => $output_name,
                      varscan_out_dir => $varscan_out_dir,
                      ref_sample_pileup => $ref_sample_pileup,
                      tumor_sample_pileup => $tumor_sample_pileup,
                      dirs => $somatic_dirs,
                      runName => $runName,
                      opt => \%opt);

        my $qsub = qsubJava(\%opt, "VARSCAN");
        if (@running_jobs) {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
        }

        push @varscan_jobs, $job_id;
    }

    my $job_id = "VS_${tumor_sample}_" . getJobId();
    my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    my $file_test = "if [ -s $ref_sample_bam -a -s $tumor_sample_bam ";
    my $snp_concat_command = catfile($opt{VCFTOOLS_PATH}, "vcf-concat");
    my $indel_concat_command = catfile($opt{VCFTOOLS_PATH}, "vcf-concat");
    my $rm_command = "rm ";

    foreach my $chr (@chrs) {
        my $snp_output = "${somatic_name}_${chr}.snp.vcf";
        my $indel_output = "${somatic_name}_${chr}.indel.vcf";
        $file_test .= "-a -s $snp_output -a -s $indel_output ";
        $snp_concat_command .= " $snp_output";
        $indel_concat_command .= " $indel_output";
        $rm_command .= "$snp_output $indel_output ";
    }

    $file_test .= "]";
    $snp_concat_command .= " > $somatic_name.snp.vcf";
    $indel_concat_command .= " > $somatic_name.indel.vcf";

    from_template("VarscanPS.sh.tt", $bash_file,
                  varscan_out_dir => $varscan_out_dir,
                  file_test => $file_test,
                  ref_sample_pileup => $ref_sample_pileup,
                  tumor_sample_pileup => $tumor_sample_pileup,
                  somatic_name => $somatic_name,
                  rm_command => $rm_command,
                  snp_concat_command => $snp_concat_command,
                  indel_concat_command => $indel_concat_command,
                  dirs => $somatic_dirs,
                  runName => $runName,
                  opt => $opt);

    my $qsub = qsubJava(\%opt, "VARSCAN");
    if (@varscan_jobs) {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @varscan_jobs) . " $bash_file";
    } else {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

sub runFreeBayes {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $somatic_dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $runName = basename($opt{OUTPUT_DIR});
    my $freebayes_out_dir = catfile($somatic_dirs->{out}, "freebayes");
    my $freebayes_tmp_dir = catfile($freebayes_out_dir, "tmp");

    if (!-d $freebayes_out_dir) {
        make_path($freebayes_out_dir) or die "Couldn't create directory $freebayes_out_dir: $!";
    }
    if (!-d $freebayes_tmp_dir) {
        make_path($freebayes_tmp_dir) or die "Couldn't create directory $freebayes_tmp_dir: $!";
    }

    my $done_file = catfile($somatic_dirs->{log}, "freebayes.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my @chrs = @{getChromosomes($opt)};
    my @freebayes_jobs;

    foreach my $chr (@chrs) {
        my $job_id = "FB_${tumor_sample}_${chr}_" . getJobId();
        my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        my $freebayes_command = catfile($opt{FREEBAYES_PATH}, "freebayes");
        $freebayes_command .= " -f $opt{GENOME} -r $chr $opt{FREEBAYES_SETTINGS} $ref_sample_bam $tumor_sample_bam > $freebayes_out_dir/$output_name.vcf";

        my $sort_uniq_filter_command = "$opt{VCFTOOLS_PATH}/vcf-sort -c -t $freebayes_tmp_dir $freebayes_out_dir/$output_name.vcf | $opt{VCFLIB_PATH}/vcfuniq > $freebayes_out_dir/$output_name.sorted_uniq.vcf";
        my $mv_command;
        if ($opt{SOMVAR_TARGETS}) {
            $sort_uniq_filter_command .= "\n\tjava -Xmx".$opt{FREEBAYES_MEM}."G -jar $opt{GATK_PATH}/GenomeAnalysisTK.jar -T SelectVariants -R $opt{GENOME} -L $opt{SOMVAR_TARGETS} -V $freebayes_out_dir/$output_name.sorted_uniq.vcf -o $freebayes_out_dir/$output_name.sorted_uniq_targetfilter.vcf\n";
            $mv_command = "mv $freebayes_out_dir/$output_name.sorted_uniq_targetfilter.vcf $freebayes_out_dir/$output_name.vcf";
        } else {
            $mv_command = "mv $freebayes_out_dir/$output_name.sorted_uniq.vcf $freebayes_out_dir/$output_name.vcf";
        }

        from_template("Freebayes.sh.tt", $bash_file,
                      chr => $chr,
                      freebayes_command => $freebayes_command,
                      sort_uniq_filter_command => $sort_uniq_filter_command,
                      mv_command => $mv_command,
                      freebayes_out_dir => $freebayes_out_dir,
                      tumor_sample_bam => $tumor_sample_bam,
                      ref_sample_bam => $ref_sample_bam,
                      dirs => $somatic_dirs,
                      runName => $runName,
                      opt => $opt);

        my $qsub = qsubJava(\%opt, "FREEBAYES");
        if (@running_jobs) {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
        }

        push @freebayes_jobs, $job_id;
    }

    my $job_id = "FB_${tumor_sample}_" . getJobId();
    my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    my $file_test = "if [ -s $ref_sample_bam -a -s $tumor_sample_bam ";
    my $concat_command = "$opt{VCFTOOLS_PATH}/vcf-concat ";
    my $rm_command = "rm -r $freebayes_tmp_dir ";
    foreach my $chr (@chrs) {
        my $snp_output = "${somatic_name}_${chr}";
        $file_test .= "-a -s ${snp_output}.vcf ";
        $concat_command .= "${snp_output}.vcf ";
        $rm_command .= "$snp_output.vcf ";
    }
    $file_test .= "]";
    $concat_command .= "> $somatic_name.vcf";

    from_template("FreebayesPostProcess.sh.tt", $bash_file,
                  file_test => $file_test,
                  concat_command => $concat_command,
                  rm_command => $rm_command,
                  ref_sample_bam => $ref_sample_bam,
                  tumor_sample_bam => $tumor_sample_bam,
                  freebayes_out_dir => $freebayes_out_dir,
                  somatic_name => $somatic_name,
                  dirs => $somatic_dirs,
                  runName => $runName,
                  opt => $opt);

    my $qsub = qsubJava(\%opt, "FREEBAYES");
    if (@freebayes_jobs) {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @freebayes_jobs) . " $bash_file";
    } else {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

sub runMutect {
    my ($tumor_sample, $somatic_name, $tumor_sample_bam, $ref_sample_bam, $running_jobs, $somatic_dirs, $opt) = @_;
    my @running_jobs = @{$running_jobs};
    my %opt = %{$opt};
    my $runName = basename($opt{OUTPUT_DIR});
    my $mutect_out_dir = catfile($somatic_dirs->{out}, "mutect");
    my $mutect_tmp_dir = catfile($mutect_out_dir, "tmp");

    if (!-d $mutect_out_dir) {
        make_path($mutect_out_dir) or die "Couldn't create directory $mutect_out_dir: $!";
    }
    if (!-d $mutect_tmp_dir) {
        make_path($mutect_tmp_dir) or die "Couldn't create directory $mutect_tmp_dir: $!";
    }

    my $done_file = catfile($somatic_dirs->{log}, "mutect.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my @chrs = @{getChromosomes($opt)};
    my @mutect_jobs;

    foreach my $chr (@chrs) {
        my $job_id = "MUT_${tumor_sample}_${chr}_" . getJobId();
        my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");
        my $output_name = "${somatic_name}_${chr}";

        my $command = "java -Xmx".$opt{MUTECT_MEM}."G -jar $opt{MUTECT_PATH}/mutect.jar -T MuTect ";
        $command .= "-R $opt{GENOME} --cosmic $opt{MUTECT_COSMIC} --dbsnp $opt{CALLING_DBSNP} --intervals $chr ";
        #if ($opt{SOMVAR_TARGETS}) { $command .= "--intervals $opt{SOMVAR_TARGETS} "; }
        $command .= "--input_file:normal $ref_sample_bam --input_file:tumor $tumor_sample_bam ";
        $command .= "--out ${output_name}.out --vcf ${output_name}_mutect.vcf";

        from_template("Mutect.sh.tt", $bash_file,
                      command => $command,
                      chr => $chr,
                      mutect_tmp_dir => $mutect_tmp_dir ,
                      tumor_sample_bam => $tumor_sample_bam,
                      ref_sample_bam => $ref_sample_bam,
                      dirs => $somatic_dirs,
                      runName => $runName,
                      opt => \%opt);

        my $qsub = qsubJava(\%opt, "MUTECT");
        if (@running_jobs) {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
        } else {
            system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
        }
        push @mutect_jobs, $job_id;
    }

    my $job_id = "MUT_${tumor_sample}_" . getJobId();
    my $bash_file = catfile($somatic_dirs->{job}, "${job_id}.sh");

    my $file_test = "if [ -s $ref_sample_bam -a -s $tumor_sample_bam ";
    my $concat_command = catfile($opt{VCFTOOLS_PATH}, "vcf-concat");
    my $filter_command = "cat ${somatic_name}_mutect.vcf | java -Xmx".$opt{MUTECT_MEM}."G -jar $opt{SNPEFF_PATH}/SnpSift.jar filter \"(na FILTER ) | (FILTER = 'PASS')\" > ${somatic_name}_mutect_passed.vcf \n";

    foreach my $chr (@chrs) {
        my $output = "${somatic_name}_${chr}_mutect.vcf";
        $file_test .= "-a -s $output ";
        $concat_command .= " $output";
    }
    $file_test .= "]";
    $concat_command .= " > ${somatic_name}_mutect.vcf";

    from_template("MutectCF.sh.tt", $bash_file,
                  mutect_tmp_dir => $mutect_tmp_dir,
                  file_test => $file_test,
                  ref_sample_bam => $ref_sample_bam,
                  tumor_sample_bam => $tumor_sample_bam,
                  somatic_name => $somatic_name,
                  concat_command => $concat_command,
                  filter_command => $filter_command,
                  mutect_out_dir => $mutect_out_dir,
                  dirs => $somatic_dirs,
                  runName => $runName,
                  opt => $opt);

    my $qsub = qsubJava(\%opt, "MUTECT");
    if (@mutect_jobs) {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id -hold_jid " . join(",", @mutect_jobs) . " $bash_file";
    } else {
        system "$qsub -o $somatic_dirs->{log} -e $somatic_dirs->{log} -N $job_id $bash_file";
    }

    return $job_id;
}

############
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
############

1;
