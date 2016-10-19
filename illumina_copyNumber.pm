package illumina_copyNumber;

use FindBin;
use lib "$FindBin::Bin";
use discipline;

use File::Basename;
use File::Spec::Functions;
use File::Path qw(make_path);

use illumina_sge qw(qsubTemplate);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadataParser qw(metadataParse);


sub runCopyNumberTools {
    my ($opt) = @_;

    my $check_cnv_jobs = [];
    if ($opt->{CNV_MODE} eq "sample_control") {
        my $metadata = metadataParse($opt->{OUTPUT_DIR});
        my $ref_sample = $metadata->{ref_sample} or die "metadata missing ref_sample";
        my $tumor_sample = $metadata->{tumor_sample} or die "metadata missing tumor_sample";

        $opt->{BAM_FILES}->{$ref_sample} or die "metadata ref_sample $ref_sample not in BAM file list";
        $opt->{BAM_FILES}->{$tumor_sample} or die "metadata tumor_sample $tumor_sample not in BAM file list";

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
    my %freec_dirs = (
        out => $out_dir,
        log => catfile($out_dir, "logs"),
        job => catfile($out_dir, "jobs"),
    );

    make_path(values %freec_dirs, { error => \my $errors });
    my $messages = join ", ", map { join ": ", each $_ } @{$errors};
    die "Couldn't create copy number output directories: $messages" if $messages;

    my @running_jobs;
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$sample}} if @{$opt->{RUNNING_JOBS}->{$sample}};
    push @running_jobs, @{$opt->{RUNNING_JOBS}->{$control}} if $control and @{$opt->{RUNNING_JOBS}->{$control}};
    my $sample_bam = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
    my $control_bam = $control ? catfile($opt->{OUTPUT_DIR}, $control, "mapping", $opt->{BAM_FILES}->{$control}) : "";

    say "\n$cnv_name \t $control_bam \t $sample_bam";

    my $done_file = catfile($freec_dirs{log}, "${cnv_name}.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my @cnv_jobs;
    if ($opt->{CNV_FREEC} eq "yes") {
        say "\n###SCHEDULING FREEC####";
        my $freec_job = runFreec($sample, $sample_bam, $control_bam, \@running_jobs, \%freec_dirs, $opt);
        push @cnv_jobs, $freec_job if $freec_job;
    }

    # check is separated from run, could have more/different CNV tools
    my $job_id = "CnvCheck_${sample}_" . getJobId();
    my $bash_file = catfile($freec_dirs{job}, "${job_id}.sh");

    from_template("CnvCheck.sh.tt", $bash_file,
                  cnv_name => $cnv_name,
                  dirs => \%freec_dirs,
                  opt => $opt);

    my $qsub = qsubTemplate($opt, "CNVCHECK");
    if (@cnv_jobs) {
        system "$qsub -o $freec_dirs{log} -e $freec_dirs{log} -N $job_id -hold_jid " . join(",", @cnv_jobs) ." $bash_file";
    } else {
        system "$qsub -o $freec_dirs{log} -e $freec_dirs{log} -N $job_id $bash_file";
    }
    push @{$check_cnv_jobs}, $job_id;
    return;
}

sub runFreec {
    my ($sample_name, $sample_bam, $control_bam, $running_jobs, $dirs, $opt) = @_;

    my $done_file = catfile($dirs->{log}, "freec.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    $dirs->{freec}{out} = catfile($dirs->{out}, "freec");
    if (!-d $dirs->{freec}{out}) {
        make_path($dirs->{freec}{out}) or die "Couldn't create directory $dirs->{freec}{out}: $!";
    }

    my @mappabilityTracks;
    @mappabilityTracks = split '\t', $opt->{FREEC_MAPPABILITY_TRACKS} if $opt->{FREEC_MAPPABILITY_TRACKS};
    my $config_file = catfile($dirs->{freec}{out}, "freec_config.txt");
    from_template("FreecConfig.tt", $config_file,
                  sample_bam => $sample_bam,
                  control_bam => $control_bam,
                  mappabilityTracks => \@mappabilityTracks,
                  dirs => $dirs,
                  opt => $opt);

    my $job_id = "Freec_${sample_name}_" . getJobId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
    my $sample_bam_name = fileparse($sample_bam);

    from_template("Freec.sh.tt", $bash_file,
                  sample_name => $sample_name,
                  sample_bam => $sample_bam,
                  control_bam => $control_bam,
                  sample_bam_name => $sample_bam_name,
                  config_file => $config_file,
                  dirs => $dirs,
                  opt => $opt);

    my $qsub = qsubTemplate($opt, "FREEC");
    if (@{$running_jobs}) {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id -hold_jid " . join(",", @{$running_jobs}) . " $bash_file";
    } else {
        system "$qsub -o $dirs->{log} -e $dirs->{log} -N $job_id $bash_file";
    }
    return $job_id;
}

1;
