package HMF::Pipeline::Finalize;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Config qw(allRunningJobs createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(fromTemplate);
use HMF::Pipeline::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run releaseLock);


sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING FINALIZE ###";

    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $extras_tar = catfile($dirs->{out}, "$opt->{RUN_NAME}_extras.tar.gz");
    my $extras_zip = catfile($dirs->{out}, "$opt->{RUN_NAME}_extras.zip");

    my $job_id = fromTemplate(
        "Finalize",
        undef,
        0,
        qsubTemplate($opt, "FINALIZE"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        done_files => $opt->{DONE_FILES},
        extras_tar => $extras_tar,
        extras_zip => $extras_zip,
        log_file => catfile($dirs->{log}, "PipelineCheck.log"),
    );

    push @{$opt->{RUNNING_JOBS}->{finalize}}, $job_id;

    HMF::Pipeline::Metadata::linkArtefact($extras_tar, "extras_tar", $opt) if $opt->{EXTRAS};
    HMF::Pipeline::Metadata::linkArtefact($extras_zip, "extras_zip", $opt) if $opt->{EXTRAS};

    return;
}

sub releaseLock {
    my ($opt) = @_;
    my $lock_file = catfile($opt->{OUTPUT_DIR}, "run.lock");
    my $dirs = createDirs($opt->{OUTPUT_DIR});
    my $job_id = fromTemplate(
        #<<< no perltidy
        "ReleaseLock",
        undef,
        0,
        qsubTemplate($opt, "FINALIZE"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        lock_file => $lock_file,
        #>>> no perltidy
    );
    return;
}

1;
