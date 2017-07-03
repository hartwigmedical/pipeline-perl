package HMF::Pipeline::DamageEstimate;

use FindBin::libs;
use discipline;

use File::Spec::Functions;

use HMF::Pipeline::Job qw(fromTemplate checkReportedDoneFile markDone);
use HMF::Pipeline::Config qw(createDirs sampleControlBamsAndJobs);
use HMF::Pipeline::Sge qw(qsubJava);

use parent qw(Exporter);
our @EXPORT_OK = qw(run);


sub run {
    my ($opt) = @_;

    my ($ref_sample, $tumor_sample, $ref_bam_path, $tumor_bam_path, $joint_name, $running_jobs) = sampleControlBamsAndJobs($opt);
    my $dirs = createDirs($opt->{OUTPUT_DIR}, damageEstimate => "damageEstimate");

    say "\n### SCHEDULING DAMAGE ESTIMATE ###";

    my $ref_job_id = fromTemplate(
        "DamageEstimate",
        undef,
        1,
        qsubJava($opt, "DAMAGE_ESTIMATE"),
        $running_jobs,
        $dirs,
        $opt,
        damage_estimate_bam_path => $ref_bam_path,
        joint_name => $joint_name
    );

    my $tumor_job_id = fromTemplate(
        "DamageEstimate",
        undef,
        1,
        qsubJava($opt, "DAMAGE_ESTIMATE"),
        $running_jobs,
        $dirs,
        $opt,
        damage_estimate_bam_path => $tumor_bam_path,
        joint_name => $joint_name
    );

    return;
}

1;
