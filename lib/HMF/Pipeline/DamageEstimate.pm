package HMF::Pipeline::DamageEstimate;

use FindBin::libs;
use discipline;
use File::Spec::Functions;
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Config qw(createDirs sampleBamAndJobs);
use HMF::Pipeline::Functions::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING DAMAGE ESTIMATE ###";

    foreach my $sample (keys %{$opt->{SAMPLES}}) {
        my $dirs = createDirs(catfile($opt->{OUTPUT_DIR}, $sample), damageEstimate => "damageEstimate");
        my ($bam_path, $running_jobs) = sampleBamAndJobs($sample, $opt);
        my $damage_estimate_job = runDamageEstimate($sample, $bam_path, $running_jobs, $dirs, $opt);

        push @{$opt->{RUNNING_JOBS}->{damageestimate}}, $damage_estimate_job;
    }
    return;
}

sub runDamageEstimate {
    my ($sample, $bam_path, $running_jobs, $dirs, $opt) = @_;
    my $out_dir_path = $dirs->{'damageEstimate'};

    say "Running damageEstimate for: $sample";
    my $job_id = fromTemplate(
        "DamageEstimate",
        $sample,
        1,
        qsubJava($opt, "DAMAGE_ESTIMATE"),
        $running_jobs,
        $dirs,
        $opt,
        sample => $sample,
        damage_estimate_bam_path => $bam_path,
        damage_estimate_out_path => $out_dir_path,
        step_name => "DamageEstimate_" . $sample,
    );

    return $job_id;
}

1;
