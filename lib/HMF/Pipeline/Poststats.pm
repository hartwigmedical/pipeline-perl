package HMF::Pipeline::Poststats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions qw(:ALL);

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING POSTSTATS ###";

    my %snp_check_dirs = map { $_ => catfile($_, "snpcheck") } keys %{$opt->{SAMPLES}};
    my $dirs = createDirs($opt->{OUTPUT_DIR}, out => "QCStats", %snp_check_dirs);
    # weird relative path requirement
    $dirs->{exoncov} = abs2rel(catfile($dirs->{out}, "exoncov"), $dirs->{tmp});

    my $done_file = catfile($dirs->{log}, "PostStats.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping";
        return;
    }

    my $sample_bams = {};
    my @running_jobs;
    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        $sample_bams->{$sample} = catfile($opt->{OUTPUT_DIR}, $sample, "mapping", $opt->{BAM_FILES}->{$sample});
        if (@{$opt->{RUNNING_JOBS}->{$sample}}) {
            push @running_jobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}});
        }
    }

    my @designs;
    @designs = split '\t', $opt->{SNPCHECK_DESIGNS} if $opt->{SNPCHECK_DESIGNS};

    my $job_id = "PostStats_" . getId();
    my $job_id_check = "PostStatsCheck_" . getId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
    my $stdout = catfile($dirs->{log}, "PostStats_$opt->{RUN_NAME}.out");
    my $stderr = catfile($dirs->{log}, "PostStats_$opt->{RUN_NAME}.err");

    writeFromTemplate(
        "PostStats.sh.tt", $bash_file,
        sample_bams => $sample_bams,
        designs => \@designs,
        job_id => $job_id,
        job_id_check => $job_id_check,
        dirs => $dirs,
        opt => $opt,
    );

    my $qsub = qsubTemplate($opt, "POSTSTATS");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    $bash_file = catfile($dirs->{job}, "${job_id_check}.sh");
    writeFromTemplate(
        "PostStatsCheck.sh.tt", $bash_file,
        designs => \@designs,
        sample_bams => $sample_bams,
        dirs => $dirs,
        opt => $opt,
    );
    system "$qsub -o $stdout -e $stderr -N $job_id_check -hold_jid $job_id $bash_file";

    $opt->{RUNNING_JOBS}->{poststats} = [$job_id_check];

    # dependent on implicit bamMetrics naming
    my $metrics_path = catfile($dirs->{out}, "$opt->{RUN_NAME}.bamMetrics.pdf");
    HMF::Pipeline::Metadata::linkArtefact($metrics_path, "qc", $opt);

    return;
}

1;
