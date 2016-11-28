package HMF::Pipeline::Poststats;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions qw(:ALL);

use HMF::Pipeline::Config qw(createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(getId);
use HMF::Pipeline::Template qw(writeFromTemplate);
use HMF::Pipeline::Metadata qw(linkArtefact linkExtraArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING POSTSTATS ###";

    my %snp_check_dirs = map { $_ => catfile($_, "snpcheck") } keys %{$opt->{SAMPLES}};
    my $dirs = createDirs($opt->{OUTPUT_DIR}, out => "QCStats", %snp_check_dirs);
    # weird relative path requirement
    my $exoncov_dir = catfile($dirs->{out}, "exoncov");
    $dirs->{exoncov} = abs2rel($exoncov_dir, $dirs->{tmp});

    my $job_id = "PostStats_" . getId();
    my $done_file = catfile($dirs->{log}, "PostStats.done");
    if (-f $done_file) {
        say "WARNING: $done_file exists, skipping $job_id";
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

    linkArtefacts(\@designs, $exoncov_dir, $dirs, $opt);
    return;
}

sub linkArtefacts {
    my ($designs, $exoncov_dir, $dirs, $opt) = @_;

    # dependent on implicit bamMetrics naming
    my $metrics_path = catfile($dirs->{out}, "$opt->{RUN_NAME}.bamMetrics.pdf");
    linkArtefact($metrics_path, "qc", $opt);

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        # dependent on implicit SNPcheck naming (easier to fix)
        foreach my $design (@{$designs}) {
            linkExtraArtefact(catfile($dirs->{$sample}, "${sample}_${design}.vcf"), $opt);
        }

        if ($opt->{EXONCALLCOV} eq "yes") {
            # dependent on implicit ExonCov naming
            linkExtraArtefact(catfile($exoncov_dir, "${sample}.html"), $opt);
            (my $coverage_name = $opt->{BAM_FILES}->{$sample}) =~ s/\.bam$/_exon_coverage.tsv/;
            linkExtraArtefact(catfile($exoncov_dir, "Exons", $coverage_name), $opt);
        }
    }
    return;
}

1;
