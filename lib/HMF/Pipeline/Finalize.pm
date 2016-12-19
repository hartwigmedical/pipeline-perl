package HMF::Pipeline::Finalize;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Config qw(allRunningJobs createDirs);
use HMF::Pipeline::Sge qw(qsubTemplate);
use HMF::Pipeline::Job qw(fromTemplate getId);
use HMF::Pipeline::Metadata;

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
        qsubTemplate($opt, "FINALIZE"),
        allRunningJobs($opt),
        $dirs,
        $opt,
        joint_name => jointName($opt),
        extras_tar => $extras_tar,
        extras_zip => $extras_zip,
        log_file => catfile($dirs->{log}, "PipelineCheck.log"),
    );

    HMF::Pipeline::Metadata::linkArtefact($extras_tar, "extras_tar", $opt) if $opt->{EXTRAS};
    HMF::Pipeline::Metadata::linkArtefact($extras_zip, "extras_zip", $opt) if $opt->{EXTRAS};

    return;
}

sub jointName {
    my ($opt) = @_;

    # do not depend on metadata if not required
    if (   $opt->{SOMATIC_VARIANTS} eq "yes"
        || ($opt->{COPY_NUMBER} eq "yes" && $opt->{CNV_MODE} eq "sample_control")
        || ($opt->{SV_CALLING} eq "yes" && $opt->{SV_MODE} eq "sample_control")) {
        my (undef, undef, $joint_name) = HMF::Pipeline::Metadata::sampleControlNames($opt);
        return $joint_name;
    }
    return "";
}

1;
