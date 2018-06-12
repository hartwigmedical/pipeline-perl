package HMF::Pipeline::HealthCheck;

use FindBin::libs;
use discipline;

use File::Find::Rule;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Sge qw(qsubJava);
use HMF::Pipeline::Functions::Config qw(allRunningJobs createDirs);

use parent qw(Exporter);
our @EXPORT_OK = qw(run
);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING HEALTH CHECKER ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $health_check_file = catfile($dirs->{log}, "$opt->{RUN_NAME}.healthcheck.json");

    my $health_check_job_id = fromTemplate(
        "HealthCheck",
        undef,
        0,
        qsubJava($opt, "HEALTHCHECK"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        run_dir => $opt->{OUTPUT_DIR},
        report_file_path => $health_check_file,
    );

    push @{$opt->{RUNNING_JOBS}->{healthcheck}}, $health_check_job_id;

    my $health_check_evaluation_job_id = fromTemplate(
        #<<< no perltidy
        "HealthCheckEvaluation",
        undef,
        0,
        qsubJava($opt, "HEALTHCHECK"),
        [ $health_check_job_id ],
        $dirs,
        $opt,
        #>>> no perltidy
    );

    push @{$opt->{RUNNING_JOBS}->{healthcheck}}, $health_check_evaluation_job_id;

    return;
}

1;
