package HMF::Pipeline::PipelineCheck;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(allRunningJobs createDirs);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Job qw(fromTemplate);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING PIPELINE CHECK ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $pipeline_check_file = "PipelineCheck.log";

    my $pipeline_check_job_id = fromTemplate(
        "PipelineCheck",
        undef,
        0,
        qsubTemplate($opt, "PIPELINE_CHECK"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        done_files => $opt->{DONE_FILES},
        log_file => catfile($dirs->{log}, $pipeline_check_file),
    );

    push @{$opt->{RUNNING_JOBS}->{pipelinecheck}}, $pipeline_check_job_id;

    return;
}

1;
