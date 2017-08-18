package HMF::Pipeline::HealthCheck;

use FindBin::libs;
use discipline;
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Config qw(allRunningJobs createDirs);

use parent qw(Exporter);
our @EXPORT_OK = qw(run
);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING HEALTH CHECKER ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});

    my $job_id = fromTemplate(
        "HealthCheck",
        undef,
        0,
        qsubJava($opt, "HEALTHCHECK"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        run_dir => $opt->{OUTPUT_DIR},
        output_path => $dirs->{log},
    );

    # TODO: run qc scripts

    push @{$opt->{RUNNING_JOBS}->{healthcheck}}, $job_id;

    return;
}

1;
