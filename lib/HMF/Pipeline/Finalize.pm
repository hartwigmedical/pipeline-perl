package HMF::Pipeline::Finalize;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(allRunningJobs createDirs);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Metadata;

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

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
        pipeline_check_file => catfile($dirs->{log}, $opt->{PIPELINE_CHECK_FILE}),
    );

    HMF::Pipeline::Metadata::linkArtefact($extras_tar, "extras_tar", $opt) if $opt->{EXTRAS};
    HMF::Pipeline::Metadata::linkArtefact($extras_zip, "extras_zip", $opt) if $opt->{EXTRAS};

    return;
}

1;
