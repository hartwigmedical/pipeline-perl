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
    my $health_check_dir = $dirs->{log};

    my $job_id = fromTemplate(
        "HealthCheck",
        undef,
        0,
        qsubJava($opt, "HEALTHCHECK"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        run_dir => $opt->{OUTPUT_DIR},
        output_path => $health_check_dir,
    );

    push @{$opt->{RUNNING_JOBS}->{healthcheck}}, $job_id;

    # TODO: run qc scripts

    my @health_check_files = File::Find::Rule->file()->name("*_health_checks.json")->in($health_check_dir);
    my $health_check_file = $health_check_files[0];
    HMF::Pipeline::Metadata::linkExtraArtefact(catfile($health_check_dir, $health_check_file), $opt);

    return;
}

1;
