package HMF::Pipeline::AmberBAFSegmentation;

use FindBin::libs;
use discipline;

use File::Basename;
use File::Spec::Functions;

use HMF::Pipeline::Functions::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Functions::Job qw(fromTemplate);
use HMF::Pipeline::Functions::Sge qw(qsubTemplate);
use HMF::Pipeline::Functions::Metadata qw(linkArtefact);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);

sub run {
    my ($opt) = @_;

    say "\n### SCHEDULING AMBER BAF SEGMENTATION ###";
    $opt->{RUNNING_JOBS}->{amber} = [];

    my $sub_dir = "amber";
    my $dirs = createDirs($opt->{OUTPUT_DIR}, amber => $sub_dir);

    my (undef, $tumor_sample, undef, undef, undef, $running_jobs) = sampleControlBamsAndJobs($opt);
    $opt->{AMBER_BAF_FILE} = "${sub_dir}/${tumor_sample}.amber.baf";
    linkArtefact($opt->{AMBER_BAF_FILE}, 'amber_baf', $opt);

    my $segmentation_job_id = fromTemplate(
        "AmberBAFSegmentation",
        undef,
        1,
        qsubTemplate($opt, "AMBER"),
        $running_jobs,
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        baf_path => $opt->{AMBER_BAF_FILE},
    );

    push @{$opt->{RUNNING_JOBS}->{amber}}, $segmentation_job_id;

    return;
}

1;
