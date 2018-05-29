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

    fromTemplate(
        "Finalize",
        undef,
        0,
        qsubTemplate($opt, "FINALIZE"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        done_files => $opt->{DONE_FILES},
        pipeline_check_file => catfile($dirs->{log}, $opt->{PIPELINE_CHECK_FILE}),
    );

    return;
}

1;
