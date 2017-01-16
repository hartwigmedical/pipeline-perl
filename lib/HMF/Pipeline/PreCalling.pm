package HMF::Pipeline::PreCalling;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(fromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    if (   ($opt->{SOMATIC_VARIANTS} eq "yes" and $opt->{SOMVAR_VARSCAN} eq "yes")
        or ($opt->{COPY_NUMBER} eq "yes" and $opt->{CNV_FREEC} eq "yes" and $opt->{FREEC_BAF} eq "yes")) {
        foreach my $sample (keys %{$opt->{SAMPLES}}) {
            my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample));
            my ($bam_path, $running_jobs) = sampleBamAndJobs($sample, $opt);
            my ($job_id, $pileup_path) = runPileup($sample, $bam_path, $running_jobs, $dirs, $opt);
            $opt->{PILEUP_FILES}->{$sample} = $pileup_path;
            push @{$opt->{RUNNING_JOBS}->{pileup}}, $job_id;
        }
    }

    return;
}

sub runPileup {
    my ($sample, $bam_path, $running_jobs, $dirs, $opt) = @_;

    say "Creating pileup for: $sample";

    (my $pileup_path = $bam_path) =~ s/\.bam$/.pileup.gz/;
    my $job_id = fromTemplate(
        "Pileup",
        $sample,
        0,
        qsubTemplate($opt, "PILEUP"),
        $running_jobs,
        $dirs,
        $opt,
        sample => $sample,
        bam_path => $bam_path,
        pileup_path => $pileup_path,
    );
    return ($job_id, $pileup_path);
}

1;
