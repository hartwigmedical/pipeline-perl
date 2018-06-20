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

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $baf_path = "${sub_dir}/${tumor_sample}.amber.baf";
    linkArtefact($baf_path, 'amber_baf', $opt);

    my $segmentation_job_id = fromTemplate(
        "AmberBAFSegmentation",
        undef,
        1,
        qsubTemplate($opt, "AMBER"),
        $running_jobs,
        $dirs,
        $opt,
        tumor_sample => $tumor_sample,
        baf_path => $baf_path,
    );

    push @{$opt->{RUNNING_JOBS}->{amber}}, $segmentation_job_id;

    return;
}

1;
