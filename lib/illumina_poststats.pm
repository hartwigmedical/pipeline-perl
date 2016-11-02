package illumina_poststats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions qw(:ALL);
use File::Path qw(make_path);

use illumina_sge qw(qsubTemplate);
use illumina_jobs qw(getJobId);
use illumina_template qw(from_template);
use illumina_metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(runPostStats);


sub runPostStats {
    my ($opt) = @_;

    say "\n### SCHEDULING POSTSTATS ###";

    my $out_dir = catfile($opt->{OUTPUT_DIR}, "QCStats");
    my $dirs = {
        out => $out_dir,
        tmp => catfile($opt->{OUTPUT_DIR}, "tmp"),
        log => catfile($opt->{OUTPUT_DIR}, "logs"),
        job => catfile($opt->{OUTPUT_DIR}, "jobs"),
    };
    map { $dirs->{$_} = catfile($dirs->{out}, $_, "snpcheck") } keys %{$opt->{SAMPLES}};

    make_path(values %{$dirs}, { error => \my $errors });
    my $messages = join ", ", map { join ": ", each $_ } @{$errors};
    die "Couldn't create output directories: $messages" if $messages;

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
            push(@running_jobs, join(",", @{$opt->{RUNNING_JOBS}->{$sample}}));
        }
    }

    my @designs;
    @designs = split '\t', $opt->{SNPCHECK_DESIGNS} if $opt->{SNPCHECK_DESIGNS};

    my $job_id = "PostStats_" . getJobId();
    my $job_id_check = "PostStatsCheck_" . getJobId();
    my $bash_file = catfile($dirs->{job}, "${job_id}.sh");
    my $stdout = catfile($dirs->{log}, "PostStats_$opt->{RUN_NAME}.out");
    my $stderr = catfile($dirs->{log}, "PostStats_$opt->{RUN_NAME}.err");

    from_template("PostStats.sh.tt", $bash_file,
                  sample_bams => $sample_bams,
                  designs => \@designs,
                  job_id => $job_id,
                  job_id_check => $job_id_check,
                  dirs => $dirs,
                  opt => $opt);

    my $qsub = qsubTemplate($opt, "POSTSTATS");
    if (@running_jobs) {
        system "$qsub -o $stdout -e $stderr -N $job_id -hold_jid " . join(",", @running_jobs) . " $bash_file";
    } else {
        system "$qsub -o $stdout -e $stderr -N $job_id $bash_file";
    }

    $bash_file = catfile($dirs->{job}, "${job_id_check}.sh");
    from_template("PostStatsCheck.sh.tt", $bash_file,
                  designs => \@designs,
                  sample_bams => $sample_bams,
                  dirs => $dirs,
                  opt => $opt);
    system "$qsub -o $stdout -e $stderr -N $job_id_check -hold_jid $job_id $bash_file";

    $opt->{RUNNING_JOBS}->{poststats} = [$job_id_check];

    # dependent on implicit bamMetrics naming
    my $metrics_path = catfile($opt->{OUTPUT_DIR}, "QCStats", "$opt->{RUN_NAME}.bamMetrics.pdf");
    illumina_metadata::linkArtefact($metrics_path, "qc", $opt);

    return;
}

1;
