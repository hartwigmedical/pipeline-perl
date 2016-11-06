package HMF::Pipeline::CopyNumber;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs addSubDir);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING COPY NUMBER TOOLS ###";

    my $check_cnv_jobs = [];
    if ($opt->{CNV_MODE} eq "sample_control") {
        my $metadata = HMF::Pipeline::Metadata::parse($opt);
        my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
        my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";

        $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};
        $opt->{BAM_FILES}->{$tumor_sample} or die "metadata tumor_sample $tumor_sample not in BAM file list: " . join ", ", keys %{$opt->{BAM_FILES}};

        my $cnv_name = "${ref_sample}_${tumor_sample}";
        runSampleCnv($tumor_sample, $ref_sample, $cnv_name, $check_cnv_jobs, $opt);
    } elsif ($opt->{CNV_MODE} eq "sample") {
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            runSampleCnv($sample, undef, $sample, $check_cnv_jobs, $opt);
        }
    }
    $opt->{RUNNING_JOBS}->{cnv} = $check_cnv_jobs;
    return;
}

sub runSampleCnv {
    my ($sample, $control, $cnv_name, $check_cnv_jobs, $opt) = @_;

    my $out_dir = catfile($opt->{OUTPUT_DIR}, "copyNumber", $cnv_name);
    my $dirs = createDirs($out_dir);

    my $running_jobs = [];
    push @{$running_jobs}, @{$opt->{RUNNING_JOBS}->{$sample}} if @{$opt->{RUNNING_JOBS}->{$sample}};
    push @{$running_jobs}, @{$opt->{RUNNING_JOBS}->{$control}} if $control and @{$opt->{RUNNING_JOBS}->{$control}};
    my $sample_bam = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
    my $control_bam = $control ? catfile($opt->{OUTPUT_DIR}, $control, "mapping", $opt->{BAM_FILES}->{$control}) : "";

    say "\n$cnv_name \t $control_bam \t $sample_bam";

    my $job_id = "CnvCheck_${sample}_" . getId();
    my $done_file = catfile($dirs->{log}, "${cnv_name}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    my @cnv_jobs;
    if ($opt->{CNV_FREEC} eq "yes") {
        say "\n###SCHEDULING FREEC####";
        my $freec_job = runFreec($sample, $sample_bam, $control_bam, $running_jobs, $dirs, $opt);
        push @cnv_jobs, $freec_job if $freec_job;
    }
    if ($opt->{CNV_QDNASEQ} eq "yes") {
        say "\n###SCHEDULING QDNASEQ####";
        my $qdnaseq_job = runQDNAseq($sample, $sample_bam, $running_jobs, $dirs, $opt);
        push @cnv_jobs, $qdnaseq_job if $qdnaseq_job;
    }

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "CnvCheck.sh.tt", $bash_file,
        cnv_name => $cnv_name,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "CNVCHECK");
    if (@cnv_jobs) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @cnv_jobs) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    push @{$check_cnv_jobs}, $job_id;
    return;
}

sub runFreec {
    my ($sample_name, $sample_bam, $control_bam, $running_jobs, $dirs, $opt) = @_;

    my $job_id = "Freec_${sample_name}_" . getId();
    my $done_file = catfile($dirs->{log}, "freec.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    $dirs->{freec}{out} = addSubDir($dirs, "freec");

    my @mappabilityTracks;
    @mappabilityTracks = split '\t', $opt->{FREEC_MAPPABILITY_TRACKS} if $opt->{FREEC_MAPPABILITY_TRACKS};
    my $config_file = catfile($dirs->{freec}{out}, "freec_config.txt");
    writeFromTemplate(
        "FreecConfig.tt", $config_file,
        sample_bam => $sample_bam,
        control_bam => $control_bam,
        mappabilityTracks => \@mappabilityTracks,
        dirs => $dirs,
        opt => $opt,
    );

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
    my $sample_bam_name = fileparse($sample_bam);

    writeFromTemplate(
        "Freec.sh.tt", $bash_file,
        sample_name => $sample_name,
        sample_bam => $sample_bam,
        control_bam => $control_bam,
        sample_bam_name => $sample_bam_name,
        config_file => $config_file,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "FREEC");
    if (@{$running_jobs}) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @{$running_jobs}) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

sub runQDNAseq {
    my ($sample_name, $sample_bam, $running_jobs, $dirs, $opt) = @_;

    my $job_id = "QDNAseq_${sample_name}_" . getId();
    my $done_file = catfile($dirs->{log}, "qdnaseq.done");
    if (-e $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
        return;
    }

    $dirs->{qdnaseq}{out} = addSubDir($dirs, "qdnaseq");

    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");

    writeFromTemplate(
        "QDNAseq.sh.tt", $bash_file,
        sample_name => $sample_name,
        sample_bam => $sample_bam,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "QDNASEQ");
    if (@{$running_jobs}) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @{$running_jobs}) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

1;
