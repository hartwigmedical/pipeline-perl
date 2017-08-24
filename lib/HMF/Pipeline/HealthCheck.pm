package HMF::Pipeline::HealthCheck;

use FindBin::libs;
use discipline;

use File::Find::Rule;
use File::Spec::Functions;

use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Sge qw(qsubJava);
use HMF::Pipeline::Config qw(allRunningJobs createDirs);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run
);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING HEALTH CHECKER ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $health_check_file = catfile($dirs->{log}, "$opt->{RUN_NAME}.healthcheck.json");

    my $job_id = fromTemplate(
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

    push @{$opt->{RUNNING_JOBS}->{healthcheck}}, $job_id;

    HMF::Pipeline::Metadata::linkExtraArtefact($health_check_file, $opt);

    return;
}

1;
